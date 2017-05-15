######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
miles.block.fast<-function(exp,year1,year2,season,DATADIR,FILESDIR)
{

#t0
t0<-proc.time()

#setting up main variables
BLOCKDIR=file.path(FILESDIR,exp,"Block",paste0(year1,"_",year2),season)
print(DATADIR)
ZDIR=file.path(DATADIR,exp)
dir.create(BLOCKDIR,recursive=T)
#outname=paste0(BLOCKDIR,"/Block_",exp,"_",year1,"_",year2,"_",season)
#outname2=paste0(BLOCKDIR,"/Events_",exp,"_",year1,"_",year2,"_",season)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#define model calendar: Gregorian, No-leap-year or 30-day?
leap_noleap=is.leapyear(years)
first_leap=years[which(leap_noleap)[1]]

#loading February leap year to check calendar in use
nomefile=paste(ZDIR,"/Z500_",exp,"_",first_leap,"02.nc",sep="")
field=ncdf.opener(nomefile,"zg","lon","lat")

#define the calendar
feb.len=length(field[1,1,])
if (feb.len==29) {kalendar="gregorian"}
if (feb.len==28) {kalendar="noleap"}
if (feb.len==30) {kalendar="30day"}
print(paste("IMPORTANT: Using",kalendar,"calendar for the time filtering"))

#setting time calendar
if (kalendar=="gregorian") {etime=power.date(season,year1,year2)}
if (kalendar=="noleap") {etime=power.date.no.leap(season,year1,year2)}
if (kalendar=="30day") {etime=power.date.30day(season,year1,year2)}

#setting up time length
ndays=0
totdays=length(etime$season)

# decleare main variables to be computed (considerable speed up!)
totrwb=totmeridional=totBI=Z500=totblocked=array(NA,dim=c(length(ics),length(ipsilon),totdays))
totTM90=array(NA,dim=c(length(ics),totdays))

# Davini et al. 2012: parameters to be set for blocking detection
fi0=30    			#lowest latitude to be analyzed
jump=15				#distance on which compute gradients
step0=round(jump/diff(ipsilon)[1])   #number of grid points to be used
central=which.min(abs(ipsilon-fi0)) 		#lowest starting latitude
north=central+step0				#lowest north latitude
south=central-step0					#lowest sourth latitude
maxsouth=central-2*step0
fiN=ipsilon[north]
fiS=ipsilon[south]	
range=round((90-fi0-jump)/diff(ipsilon)[1])  #escursion to the north for computing blocking

# TM90: parametres for blocking detection (beta)
tm90_fi0=60 #central_lat
tm90_fiN=tm90_fi0+20; tm90_fiS=tm90_fi0-20 #south and north lath
tm90_central=whicher(ipsilon,tm90_fi0)
tm90_south=whicher(ipsilon,tm90_fiS)
tm90_north=whicher(ipsilon,tm90_fiN)
tm90_range=c(-2:2)  #5 degrees to the north, 5 to the south



print("--------------------------------------------------")
print(c("distance for gradients:",step0*diff(ics)[1]))
print(paste("range of latitudes ",fi0,"-",90-step0*diff(ics)[1]," N",sep=""))

print("--------------------------------------------------")

##########################################################
#--------------Istantaneous Blocking---------------------#
##########################################################

# main cycles on years and months
for (yy in year1:year2)
{
for (mm in timeseason)
		{
		mm=sprintf("%02d",mm)
                nomefile=paste0(ZDIR,"/Z500_",exp,"_",yy,mm,".nc")
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
	
	
		#introduce matrix for blocking computation
		BI=rwb=meridional=array(NA,dim=c(length(ics),length(ipsilon),days))
		blocked=array(0,dim=c(length(ics),length(ipsilon),days))
		tm90_blocked=array(0,dim=c(length(ics),days))

		#----COMPUTING BLOCKING INDICES-----
		for (t in 1:days)
		{

		#TM90: beta
		tm90_ghgn=(field[,tm90_north+tm90_range,t]-field[,tm90_central+tm90_range,t])/(tm90_fiN-tm90_fi0)
                tm90_ghgs=(field[,tm90_central+tm90_range,t]-field[,tm90_south+tm90_range,t])/(tm90_fi0-tm90_fiS)
                tm90_check=(tm90_ghgs>0 & tm90_ghgn<(-10)) # TM90 conditions
                tm90_check[tm90_check==T]=1; tm90_check[tm90_check==F]=0
		tm90_blocked[,t]=apply(tm90_check,c(1),max,na.rm=T)
		
		
		#multidim extension
                new_field=rbind(field[,,t],field[,,t],field[,,t])

                for (delta in 0:range)  # computing blocking for different latitudes
                {
                ghgn=(field[,north+delta,t]-field[,central+delta,t])/(fiN-fi0)
                ghgs=(field[,central+delta,t]-field[,south+delta,t])/(fi0-fiS)
                gh2gs=(field[,south+delta,t]-field[,maxsouth+delta,t])/(fi0-fiS)
                check1=which(ghgs>0 & ghgn<(-10))
                #check1=which(ghgs>0 & ghgn<(-10) & gh2gs<(-5))  #supplementary condition

                if (length(check1)>0)
                {
                # 1-MATRIX FOR INSTANTANEOUS BLOCKING
                blocked[check1,central+delta,t]=1

                # 2-PART ON COMPUTATION OF ROSSBY WAVEBREAKING
                r=check1+length(ics)
                rwb_jump=jump/2
                steprwb=round(rwb_jump/(ics[20]-ics[19]))
                rwb_west=new_field[(r-steprwb),south+delta+steprwb]
                rwb_east=new_field[(r+steprwb),south+delta+steprwb]
                fullgh=(rwb_west-rwb_east)

                rwb[check1[fullgh<0],central+delta,t]=(-10)     # gradient decreasing: cyclonic RWB     
                rwb[check1[fullgh>0],central+delta,t]=10        # gradient increasing: anticyclonic RWB                 

                # 4-part about adapted version of blocking intensity by Wiedenmann et al. (2002)
                step=round(60/(ics[2]-ics[1]))
                ii=check1+length(ics)
                zu=zd=NULL
                for (ll in ii)
                        {
                        zu=c(zu,min(new_field[(ll-step):ll,central+delta]))
                        zd=c(zd,min(new_field[ll:(ll+step),central+delta]))
                        }
                mz=field[check1,central+delta,t]
                rc=0.5*((zu+mz)/2+(zd+mz)/2)
                BI[check1,central+delta,t]=100*(mz/rc-1)

                #5 - part about meridional gradient index
                meridional[check1,central+delta,t]=ghgs[check1]


                }}}}

		# prepare data for arrays
		daystowrite=(ndays+1):(ndays+days)
		ndays=ndays+days
	
		# update arrays
		totblocked[,,daystowrite]=blocked
		Z500[,,daystowrite]=field
		totBI[,,daystowrite]=BI
		totrwb[,,daystowrite]=rwb
		totmeridional[,,daystowrite]=meridional
		totTM90[,daystowrite]=tm90_blocked
		
                print(paste("Total # of days:",ndays))
                print("-------------------------")

}}

##########################################################
#--------------------Mean Values-------------------------#
##########################################################

# beta for TM90
TM90=apply(totTM90,1,mean)

#compute mean values (use rowMeans that is faster when there are no NA values)
frequency=rowMeans(totblocked,dims=2)*100               #frequency of Instantaneous Blocking days
Z500mean=rowMeans(Z500,dims=2)                          #Z500 mean value
BI=apply(totBI,c(1,2),mean,na.rm=T)                     #Blocking Intensity Index as Wiedenmann et al. (2002)
MGI=apply(totmeridional,c(1,2),mean,na.rm=T)   		 #Value of meridional gradient inversion

#anticyclonic and cyclonic averages RWB
CN=apply(totrwb,c(1,2),function(x) sum(x[x==(-10)],na.rm=T))/(ndays)*(-10)
ACN=apply(totrwb,c(1,2),function(x) sum(x[x==(10)],na.rm=T))/(ndays)*(10)

#saving
#print("Saving Instantaneous Blocking...")
#save(ics,ipsilon,totblocked,frequency,Z500,Z500mean,totmeridional,totBI,totrwb,ACN,CN,MGI,BI,file=outname)

t1=proc.time()-t0
print(t1)

##########################################################
#--------------------Time filtering----------------------#
##########################################################

#spatial filtering on fixed longitude distance
spatial=longitude.filter(ics,ipsilon,totblocked)
#CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

#large scale extension on 10x5 box
large=largescale.extension.if(ics,ipsilon,spatial)
#large=largescale.extension2(ics,ipsilon,spatial)
#LARGE=apply(large,c(1,2),sum,na.rm=T)/ndays*100

#5 days persistence filter
block=blocking.persistence(large,persistence=5,time.array=etime)

#saving
#print("Saving Blocking Events...")
#save(ics,ipsilon,block,file=outname2)

tf=proc.time()-t1
print(tf)


##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
savefile1=paste0(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc")
savefile2=paste0(BLOCKDIR,"/BlockFull_",exp,"_",year1,"_",year2,"_",season,".nc")

#which fieds to plot/save
fieldlist=c("InstBlock","Z500","MGI","BI","CN","ACN","BlockEvents","DurationEvents","NumberEvents")
full_fieldlist=c("InstBlock","Z500","MGI","BI","CN","ACN","BlockEvents")

for (var in fieldlist)
{
        #name of the var
        if (var=="InstBlock")
                {longvar="Instantaneous Blocking frequency"; unit="%"; field=frequency; full_field=totblocked}
        if (var=="Z500")
                {longvar="Geopotential Height"; unit="m"; field=Z500mean; full_field=Z500}
        if (var=="BI")
                {longvar="BI index"; unit=""; field=BI; full_field=totBI}
        if (var=="MGI")
                {longvar="MGI index"; unit=""; field=MGI; full_field=totmeridional}
        if (var=="ACN")
                {longvar="Anticyclonic RWB frequency"; unit="%"; field=ACN; full_field=totrwb/10; full_field[full_field==(-1)]=NA}
        if (var=="CN")
                {longvar="Cyclonic RWB frequency"; unit="%"; field=CN; full_field=totrwb/10; full_field[full_field==(1)]=NA}
        if (var=="BlockEvents")
                {longvar="Blocking Events frequency"; unit="%"; field=block$percentage; full_field=block$track}
        if (var=="DurationEvents")
                {longvar="Blocking Events duration"; unit="days"; field=block$duration}
        if (var=="NumberEvents")
                {longvar="Blocking Events number"; unit=""; field=block$nevents}

	#fix eventual NaN	
	field[is.nan(field)]=NA

        # dimensions definition
        TIME=paste("days since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
        LEVEL=50000
	fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
        x <- ncdim_def( "Lon", "degrees", ics)
        y <- ncdim_def( "Lat", "degrees", ipsilon)
        z <- ncdim_def( "Lev", "Pa", LEVEL)
        t1 <- ncdim_def( "Time", TIME, 1,unlim=T)
	t2 <- ncdim_def( "Time", TIME, fulltime,unlim=T)
	

        #variable definitions
        var_ncdf=ncvar_def(var,unit,list(x,y,z,t=t1),-999,longname=longvar,prec="single",compression=1)
	full_var_ncdf=ncvar_def(var,unit,list(x,y,z,t=t2),-999,longname=longvar,prec="single",compression=1)
	
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
for (var in fieldlist)
{
        # put variables into the ncdf file
        ncvar_put(ncfile1, fieldlist[which(var==fieldlist)], get(paste0("field",var)), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
}
nc_close(ncfile1)

#Fullfield Netcdf file creation
print(savefile2)
namelist2=paste0("full_var",full_fieldlist)
nclist2 <- mget(namelist2)
ncfile2 <- nc_create(savefile2,nclist2)
for (var in full_fieldlist)
{
        # put variables into the ncdf file
        ncvar_put(ncfile2, full_fieldlist[which(var==full_fieldlist)], get(paste0("full_field",var)), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
}
nc_close(ncfile2)

}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","season","DATADIR","FILESDIR","PROGDIR")
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
        miles.block.fast(exp,year1,year2,season,DATADIR,FILESDIR)
    }
}
