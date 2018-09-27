######################################################
#-----Meandering routines computation for MiLES------#
#--------G. Di Capua & P. Davini (Apr 2018)----------#
######################################################
miles.meandering<-function(dataset,expid,ens,year1,year2,season,z500filename,FILESDIR,doforce) {

#t0
t0<-proc.time()

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#define folders using file.builder function (takes care of ensembles)
savefile1=file.builder(FILESDIR,"MI","MIClim",dataset,expid,ens,year1,year2,season)
savefile2=file.builder(FILESDIR,"MI","MIFull",dataset,expid,ens,year1,year2,season)

#check if data is already there to avoid re-run
if (file.exists(savefile1)) {
       print("Actually requested Meandering Index data is already there!")
       if (doforce=="true") {
               print("Running with doforce=true... re-run!")
       } else  {
       print("Skipping... activate doforce=true if you want to re-run it"); q() 
       }
}

fieldlist=ncdf.opener.universal(z500filename,namevar="zg",tmonths=timeseason,tyears=years,rotate="full")
print(str(fieldlist))


#extract calendar and time unit from the original file
tcal=attributes(fieldlist$time)$cal
tunit=attributes(fieldlist$time)$units

#time array to simplify time filtering
etime=power.date.new(fieldlist$time)
totdays=length(fieldlist$time)

#declare variable
Z500=fieldlist$field

# list of isohypses on which evaluate the MI
isolvls=seq(4900,6200,5)

# reference latitude (60N following Di Capua et al., 2016)
ref_lat=60

#Running the real code
MI_lat=MI_value=1:totdays*NA

#library(rbenchmark)
#print(benchmark(sapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)}), 
#		vapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)},list(1,2))))


for (t in 1:totdays) {
	progression.bar(t,totdays,each=20)
	MI_list=sapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,t],x,ref_lat,verbose=F)})
	#MI_list=vapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)},list(1,2))
        MI_value[t]=max(unlist(MI_list[1,]))
        MI_lat[t]=unlist(MI_list[2,])[which.max(unlist(MI_list[1,]))]
}


tf=proc.time()-t0
print(tf)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")

#which fieds to plot/save
fieldlist=c("MI","MI_lat")
full_fieldlist=c("MI","MI_lat")

# dimensions definition
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
TIME=paste(tunit," since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
x <- ncdim_def( "lon", "degrees_east", 0, longname="longitude")
t1 <- ncdim_def( "time", TIME, 0, unlim=T, calendar=tcal, longname="time")
t2 <- ncdim_def( "time", TIME, fulltime,unlim=T, calendar=tcal, longname="time")


for (var in fieldlist)
{
        #name of the var
        if (var=="MI")
                {longvar="Meandering Index"; unit=""; field=mean(MI_value); full_field=MI_value}
        if (var=="MI_lat")
                {longvar="Meandering Index Latitude"; unit="deg"; field=mean(MI_lat); full_field=MI_lat}

        #variable definitions
        var_ncdf=ncvar_def(var,unit,list(x,t=t1),-999,longname=longvar,prec="single",compression=1)
        full_var_ncdf=ncvar_def(var,unit,list(x,t=t2),-999,longname=longvar,prec="single",compression=1)

        assign(paste0("var",var),var_ncdf)
        assign(paste0("full_var",var),full_var_ncdf)
        assign(paste0("field",var),field)
        assign(paste0("full_field",var),full_field)
}

#Climatologies Netcdf file creation
print(savefile1)
namelist1=paste0("var",fieldlist)
nclist1 <- mget(namelist1)
ncfile1 <- nc_create(savefile1,nclist1)
for (var in fieldlist) {
        # put variables into the ncdf file
        ncvar_put(ncfile1, fieldlist[which(var==fieldlist)], get(paste0("field",var)), start = c(1, 1),  count = c(-1,-1))
}
nc_close(ncfile1)

#Fullfield Netcdf file creation
print(savefile2)
namelist2=paste0("full_var",full_fieldlist)
nclist2 <- mget(namelist2)
ncfile2 <- nc_create(savefile2,nclist2)
for (var in full_fieldlist) {
        # put variables into the ncdf file
        ncvar_put(ncfile2, full_fieldlist[which(var==full_fieldlist)], get(paste0("full_field",var)), start = c(1, 1),  count = c(-1,-1))

}
nc_close(ncfile2)

}

#blank lines
cat("\n\n\n")
# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("dataset","expid","ens","year1","year2","season","z500filename","FILESDIR","PROGDIR","doforce")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
        print(name_args)
    } else {
# when the number of arguments is ok run the function()
        for (k in 1:req_args) {assign(name_args[k],args[k])}
        source(file.path(PROGDIR,"script/basis_functions.R"))
        source(file.path(PROGDIR,"script/meandering_functions.R"))
        miles.meandering(dataset,expid,ens,year1,year2,season,z500filename,FILESDIR,doforce)
    }
}

