######################################################
#------Regimes routines computation for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#get environmental variables
#PROGDIR<-Sys.getenv(c("PROGDIR"))
#ZDIR<-Sys.getenv(c("ZDIR"))
#BLOCKDIR<-Sys.getenv(c("BLOCKDIR"))

PROGDIR="/home/paolo/MiLES"
ZDIR="/work/users/paolo/scratch/miles/Z500/ERAINTERIM"

#t0
t0<-proc.time()

#read command line
#args <- commandArgs(TRUE)
#exp=args[1]
#year1=args[2]
#year2=args[3]
#season=args[4]
exp="ERAINTERIM"
year1=1979
year2=2014
season="DJF"

#setting up main variables
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))
ZDIR=paste(ZDIR,"/",sep="")
#BLOCKDIR=paste(BLOCKDIR,"/",year1,"_",year2,"/",season,"/",sep="")
#dir.create(BLOCKDIR,recursive=T)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#define model calendar: Gregorian, No-leap-year or 30-day?
leap_noleap=is.leapyear(years)
first_leap=years[which(leap_noleap)[1]]

#loading February leap year to check calendar in use
nomefile=paste(ZDIR,"Z500_",exp,"_",first_leap,"02.nc",sep="")
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
                nomefile=paste(ZDIR,"Z500_",exp,"_",yy,mm,".nc",sep="")
                print(c(season,exp,yy,mm))

                ##check existance of the files
                if (file.exists(nomefile)==FALSE)    # if files do not exist add empty array for all fields
                {
                print("last file does not exist")
                print("File missing: aborted")
                break
                }

                if (file.exists(nomefile)==TRUE)
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

Z500mean=rowMeans(Z500,dims=2) 
Z500anom=sweep(Z500,c(1,2),Z500mean,"-")

xlim=c(-90,40)
ylim=c(20,85)
seeds=4
weather_regimes=regimes(ics,ipsilon,Z500anom,ncluster=seeds,ntime=1000,neof=10,xlim,ylim)

kk=order(weather_regimes$frequencies,decreasing=T)
frequencies=weather_regimes$frequencies[kk]
pattern=weather_regimes$regimes[,,kk]
timeseries=weather_regimes$cluster+10
for (ss in 1:seeds)
        {timeseries[timeseries==(ss+10)]=which(kk==ss)}



t1=proc.time()-t0
print(t1)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
BLOCKDIR="/work/users/paolo/scratch/"
savefile1=paste(BLOCKDIR,"/RegimesPattern_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
savefile2=paste(BLOCKDIR,"/RegimesTimeseries_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
savefile3=paste(BLOCKDIR,"/RegimesFrequencies_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")

# dimensions definition
TIME=paste("days since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
x <- ncdim_def( "Lon", "degrees", ics)
y <- ncdim_def( "Lat", "degrees", ipsilon)
z <- ncdim_def( "Lev", "Pa", LEVEL)
t0 <- ncdim_def( "Time", TIME, 1:4,unlim=T)
t1 <- ncdim_def( "Time", TIME, 1,unlim=T)
t2 <- ncdim_def( "Time", TIME, fulltime,unlim=T)

var="Regimes"; unit="m"; longvar="North Atlantic Regimes"
pattern_ncdf=ncvar_def("Regimes",unit,list(x,y,z,t=t0),-999,longname=longvar,prec="single",compression=1)

ncfile1 <- nc_create(savefile1,pattern_ncdf)
ncvar_put(ncfile1, var, pattern, start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
nc_close(ncfile1)

#map("world",regions=".",interior=F,exact=F,boundary=T,add=T)

