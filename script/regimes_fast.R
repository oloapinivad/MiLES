######################################################
#------Regimes routines computation for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles.regimes.fast<-function(exp,ens,year1,year2,season,z500filename,FILESDIR,nclusters=nclusters,doforce)
{

#t0
t0<-proc.time()

#region boundaries for North Atlantic 
if (nclusters!=4 | season!="DJF") {
	stop("Beta version: unsupported season and/or number of clusters")
}

#test function to smooth seasonal cycle: it does not work fine yet, keep it false
smoothing=F
xlim=c(-80,40)
ylim=c(30,87.5)

#define file where save data

savefile1=file.builder(FILESDIR,"Regimes","RegimesPattern",exp,ens,year1,year2,season)
#check if data is already there to avoid re-run
if (file.exists(savefile1)) {
        print("Actually requested weather regimes data is already there!")
	if (doforce=="true") {
                print("Running with doforce=true... re-run!")
        } else  {
	        print("Skipping... activate doforce=true if you want to re-run it"); q()
        }
}


#setting up main variables
#REGIMESDIR=file.path(FILESDIR,exp,"Regimes",paste0(year1,"_",year2),season)
#dir.create(REGIMESDIR,recursive=T)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#if (smoothing) {
#	timeseason0=timeseason
#	if (season=="DJF") {
#		timeseason=sort(c(timeseason,11,3)) 
#	} else {
#	timeseason=sort(c(timeseason[1]-1,timeseason,timeseason[length(timeseason)]+1))
#	}
#}

#new file opening
nomefile=z500filename
fieldlist=ncdf.opener.time(nomefile,"zg",tmonths=timeseason,tyears=years,rotate="full")

#extract calendar and time unit from the original file
tcal=attributes(fieldlist$time)$cal
tunit=attributes(fieldlist$time)$units

#time array
etime=power.date.new(fieldlist$time)

#declare variable
Z500=fieldlist$field

print("Compute anomalies based on daily mean")
#old clean script with apply and ave, slow
#Z500cycle=apply(Z500,c(1,2),ave,etime$month,etime$day)
#if (!smoothing) {
#    Z500anom=Z500-aperm(Z500cycle,c(2,3,1)) 
#}

#beta function for daily anomalies, use array predeclaration and rowMeans (40 times faster!)
daily.anom.mean.beta<-function(ics,ipsilon,field,etime) {
    condition=paste(etime$day,etime$month)
    daily=array(NA,dim=c(length(ics),length(ipsilon),length(unique(condition))))
    anom=field*NA
    for (t in unique(condition)) {
        daily[,,which(t==unique(condition))]=rowMeans(field[,,t==condition],dims=2)
        anom[,,which(t==condition)]=sweep(field[,,which(t==condition)],c(1,2),daily[,,which(t==unique(condition))],"-")
    }
    return(anom)
}

Z500anom=daily.anom.mean.beta(ics,ipsilon,Z500,etime)

#if (smoothing) {
#	print("running mean")
#	rundays=5
#	runZ500cycle=apply(Z500cycle,c(2,3),filter,rep(1/rundays,rundays),sides=2)
#	Z500anom0=Z500-aperm(runZ500cycle,c(2,3,1))
#	whichdays=which(as.numeric(format(datas,"%m")) %in% timeseason0)
#	Z500anom=Z500anom0[,,whichdays]
#	etime=power.date.new(datas[whichdays])
#	}

#compute weather regimes
weather_regimes=regimes(ics,ipsilon,Z500anom,ncluster=nclusters,ntime=1000,neof=4,xlim,ylim,alg="Hartigan-Wong")

t1=proc.time()-t0
print(t1)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")

#--deprecated--# 
#print(fulltime[2])
#temporary check for seconds/days TO BE FIXED
#if (fulltime[2]==1) {tunit="days"}
#if (fulltime[2]==86400) {tunit="seconds"}
#--------------#

# dimensions definition
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
TIME=paste(tunit," since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
x <- ncdim_def( "Lon", "degrees_east", ics, longname="Longitude")
y <- ncdim_def( "Lat", "degrees_north", ipsilon, longname="Latitude")
z <- ncdim_def( "Lev", "Pa", LEVEL, longname="Pressure")
t <- ncdim_def( "Time", TIME, fulltime,calendar=tcal, longname="Time", unlim=T)

# extra dimensions definition
x0 <- ncdim_def( "lon", "degrees_east", 0, longname="Longitude")
y0 <- ncdim_def( "lat", "degrees_north", 0, longname="Longitude")
cl <- ncdim_def( "Lev", "cluster index", 1:nclusters, longname="Pressure")

#var definition
unit="m"; longvar="Weather Regimes Pattern"
pattern_ncdf=ncvar_def("Regimes",unit,list(x,y,cl),-999,longname=longvar,prec="single",compression=1)

unit=paste0("0-",nclusters); longvar="Weather Regimes Cluster Index"
cluster_ncdf=ncvar_def("Indices",unit,list(t),-999,longname=longvar,prec="single",compression=1)

unit="%"; longvar="Weather Regimes Frequencies"
frequencies_ncdf=ncvar_def("Frequencies",unit,list(cl),-999,longname=longvar,prec="single",compression=1)

#saving file
ncfile1 <- nc_create(savefile1,list(pattern_ncdf,cluster_ncdf,frequencies_ncdf))
ncvar_put(ncfile1, "Regimes", weather_regimes$regimes, start = c(1, 1, 1),  count = c(-1,-1,-1))
ncvar_put(ncfile1, "Indices", weather_regimes$cluster, start = c(1),  count = c(-1))
ncvar_put(ncfile1, "Frequencies", weather_regimes$frequencies, start = c(1),  count = c(-1))
nc_close(ncfile1)

}

#blank line
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","ens","year1","year2","season","z500filename","FILESDIR","PROGDIR","nclusters","doforce")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
        print(name_args)
    } else {
# when the number of arguments is ok run the function()
        for (k in 1:req_args) {assign(name_args[k],args[k])}
        source(paste0(PROGDIR,"/script/basis_functions.R"))
        miles.regimes.fast(exp,ens,year1,year2,season,z500filename,FILESDIR,nclusters,doforce)
    }
}


