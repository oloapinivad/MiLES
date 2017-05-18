######################################################
#------Regimes routines computation for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles.regimes.fast<-function(exp,year1,yera2,season,DATADIR,FILESDIR,nclusters=nclusters)
{

#t0
t0<-proc.time()

#region boundaries for North Atlantic 
NorthAtlantic=T
xlim=c(-80,40)
ylim=c(30,87.5)

#setting up main variables
REGIMESDIR=file.path(FILESDIR,exp,"Regimes",paste0(year1,"_",year2),season)
dir.create(REGIMESDIR,recursive=T)
ZDIR=file.path(DATADIR,exp)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#define model calendar: Gregorian, No-leap-year or 30-day?
leap_noleap=is.leapyear(years)
first_leap=years[which(leap_noleap)[1]]

#loading February leap year to check calendar in use
nomefile=paste0(ZDIR,"/Z500_",exp,"_",first_leap,"02.nc")
field=ncdf.opener(nomefile,"zg","lon","lat")

#define the calendar
feb.len=length(field[1,1,])
if (feb.len==29) {kalendar="gregorian"}
if (feb.len==28) {kalendar="noleap"}
if (feb.len==30) {kalendar="30day"}
print(paste("IMPORTANT: Using",kalendar,"calendar"))

#setting time calendar
if (kalendar=="gregorian") {etime=power.date(season,year1,year2)}
if (kalendar=="noleap") {etime=power.date.no.leap(season,year1,year2)}
if (kalendar=="30day") {etime=power.date.30day(season,year1,year2)}

#setting up time length
totdays=length(etime$season)
ndays=0

#declare variable
Z500=array(NA,dim=c(length(ics),length(ipsilon),totdays))

# main cycles on years and months
for (yy in year1:year2)
{
for (mm in timeseason)
                {
                mm=sprintf("%02d",mm)
                nomefile=paste0(ZDIR,"/Z500_",exp,"_",yy,mm,".nc")
                print(c(season,exp,yy,mm))

                ##check existance of the files
                if (!file.exists(nomefile))    # if files do not exist add empty array for all fields
                {
                print("last file does not exist")
                stop("File missing: aborted")
                }

                if (file.exists(nomefile))
                {

                # routine for reading the file
                field=ncdf.opener(nomefile,"zg","lon","lat")
		
                days=length(field[1,1,])
		daystowrite=(ndays+1):(ndays+days)
                ndays=ndays+days
		Z500[,,daystowrite]=field
		}
		print(paste("Total # of days:",ndays))
                print("-------------------------")
		}
}

#to adapt to Irene's script
#k0=60:(60+2617)
#Z500=Z500[,,k0]

print("Compute anomalies based on daily mean")
#Z500cycle=apply(Z500[,,k0],c(1,2),ave,etime$month[k0],etime$day[k0])
Z500cycle=apply(Z500,c(1,2),ave,etime$month,etime$day)
Z500anom=Z500-aperm(Z500cycle,c(2,3,1))

#print("running mean")
#runZ500cycle=apply(Z500cycle,c(2,3),filter,rep(1/5,5),sides=2)
#runZ500cycle[1:2,,]=runZ500cycle[3,,]
#runZ500cycle[2617:2618,,]=runZ500cycle[2616,,]
#Z500anom=Z500-aperm(runZ500cycle,c(2,3,1))

#Z500anom4=ncdf.opener("/home/paolo/scratch/miles/anomaly_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_0.nc",rotate="no")
#Z500anom2=ncdf.opener("/home/paolo/scratch/miles/final.nc")
#Z500anom3=ncdf.opener("/home/paolo/scratch/miles/final_DJF.nc")


print(str(Z500anom))
weather_regimes=regimes(ics,ipsilon,Z500anom,ncluster=nclusters,ntime=1000,neof=4,xlim,ylim,alg="Hartigan-Wong")
#weather_regimes4=regimes(slon,slat,Z500anom4,ncluster=nclusters,ntime=1000,neof=4,xlim,ylim,alg="Lloyd")
#weather_regimes5=regimes(slon,slat,Z500anom4,ncluster=nclusters,ntime=1000,neof=4,xlim,ylim,alg="Hartigan-Wong")


#irene=read.table("/home/paolo/devel4miles/indclORD_4clus_zg500_day_ERAInterim_obs_144x73_1ens_DJF_EAT_4pcs.txt")$V1
#ifreq=table(irene)/length(irene)*100
#print(ifreq)

t1=proc.time()-t0
print(t1)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
savefile1=paste(REGIMESDIR,"/RegimesPattern_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")

# dimensions definition
TIME=paste("days since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
x <- ncdim_def( "Lon", "degrees", ics)
x0 <- ncdim_def( "Lon0", "degrees", 0)
y <- ncdim_def( "Lat", "degrees", ipsilon)
y0 <- ncdim_def( "Lat0", "degrees", 0)
z <- ncdim_def( "Lev", "Pa", LEVEL)
cl <- ncdim_def( "Time0", TIME, 1:nclusters)
t <- ncdim_def( "Time", TIME, fulltime,unlim=T)

unit="m"; longvar="Weather Regimes Pattern"
pattern_ncdf=ncvar_def("Regimes",unit,list(x,y,z,cl),-999,longname=longvar,prec="single",compression=1)
unit=paste0("0-",nclusters); longvar="Weather Regimes Cluster Undex"
cluster_ncdf=ncvar_def("Indices",unit,list(x0,y0,z,t),-999,longname=longvar,prec="single",compression=1)
unit="%"; longvar="Weather Regimes Frequencies"
frequencies_ncdf=ncvar_def("Frequencies",unit,list(cl),-999,longname=longvar,prec="single",compression=1)

ncfile1 <- nc_create(savefile1,list(pattern_ncdf,cluster_ncdf,frequencies_ncdf))
ncvar_put(ncfile1, "Regimes", weather_regimes$regimes, start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
ncvar_put(ncfile1, "Indices", weather_regimes$cluster, start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
ncvar_put(ncfile1, "Frequencies", weather_regimes$frequencies, start = c(1),  count = c(-1))
nc_close(ncfile1)





}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","season","DATADIR","FILESDIR","PROGDIR","nclusters")
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
        miles.regimes.fast(exp,year1,year2,season,DATADIR,FILESDIR,nclusters)
    }
}


