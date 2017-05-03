######################################################
#-----Blocking routines computation for ECE3---------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################

#get environmental variables
INDIR<-Sys.getenv(c("INDIR"))
PROGDIR<-Sys.getenv(c("PROGDIR"))
ZDIR<-Sys.getenv(c("ZDIR"))
BLOCKDIR<-Sys.getenv(c("BLOCKDIR"))
#R_LIBLOC<-Sys.getenv(c("R_LIBLOC"))
R_LIBLOC=.libPaths()[1]

#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]

#loadin packages
library("ncdf4",lib.loc=R_LIBLOC)

#setting up main variables
years=year1:year2
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))
timeseason=season2timeseason(season)
ZDIR=paste(ZDIR,"/",sep="")
BLOCKDIR=paste(BLOCKDIR,"/",year1,"_",year2,"/",season,"/",sep="")
dir.create(BLOCKDIR,recursive=T)
outname=paste(BLOCKDIR,"/Block_",exp,"_",year1,"_",year2,"_",season,sep="")

#loading first file to prepare the matricies
nomefile=paste(ZDIR,"Z500_",exp,"_",year1,"01.nc",sep="")
field=ncdf.opener(nomefile,"zg","lon","lat")
if (min(field)<10000) {gravity=1} else {gravity=9.8065}

ndays=0

# decleare main variables to be computed
totrwb=totmeridional=totBI=Z500=totblocked=array(NA,dim=c(length(ics),length(ipsilon)))

# parameters to be set for blocking detection
fi0=30    			#lowest latitude to be analyzed
jump=15				#distance on which compute gradients
step0=round(jump/(ipsilon[2]-ipsilon[1]))   #number of grid points to be used
central=which.min(abs(ipsilon-fi0)) 		#lowest starting latitude
north=central+step0				#lowest north latitude
south=central-step0					#lowest sourth latitude
maxsouth=central-2*step0
fiN=ipsilon[north]
fiS=ipsilon[south]	
range=round((90-fi0-jump)/(ipsilon[20]-ipsilon[19]))  #escursion to the north for computing blocking

print("--------------------------------------------------")
print(c("distance for gradients:",step0*2.5))
print(paste("range of latitudes ",fi0,"-",90-step0*2.5," N",sep=""))

print("--------------------------------------------------")


# main cycles on years and months
for (yy in years)
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
		print("File missing: future procedure of tracking aborted")
		break
		}
		 
		if (file.exists(nomefile)==TRUE) 
		{

		# routine for reading the file
		field=ncdf.opener(nomefile,"zg","lon","lat")
		days=length(field[1,1,])
	
		#introduce matrix for blocking computation
		BI=rwb=meridional=array(NA,dim=c(length(ics),length(ipsilon),length(field[1,1,])))
		blocked=array(0,dim=c(length(ics),length(ipsilon),length(field[1,1,])))

		#----COMPUTING BLOCKING INDICES-----
		for (t in 1:(length(field[1,1,])))
		{
		
		field[,,t]=field[,,t]/gravity
		
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
                rwb_jump=7.5
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

		# concatenate matrix
		totblocked=c(totblocked,blocked)
		Z500=c(Z500,field)
		totBI=c(totBI,BI)
		totmeridional=c(totmeridional,meridional)
		totrwb=c(totrwb,rwb)	
		
		ndays=ndays+days
                print(paste("Total # of days:",ndays))
                print("-------------------------")

}}


#reconstruct the matrix in 3D form (lon,lat,time)

totblocked=array(totblocked,dim=c(length(ics),length(ipsilon),length(totblocked)/(length(ics)*length(ipsilon))))
totmeridional=array(totmeridional,dim=c(length(ics),length(ipsilon),length(totmeridional)/(length(ics)*length(ipsilon))))
Z500=array(Z500,dim=c(length(ics),length(ipsilon),length(Z500)/(length(ics)*length(ipsilon))))
totBI=array(totBI,dim=c(length(ics),length(ipsilon),length(totBI)/(length(ics)*length(ipsilon))))
totrwb=array(totrwb,dim=c(length(ics),length(ipsilon),length(totrwb)/(length(ics)*length(ipsilon))))

#compute mean values
frequency=apply(totblocked,c(1,2),mean,na.rm=T)*100       #frequency of Instantaneous Blocking days
BI=apply(totBI,c(1,2),mean,na.rm=T)                     #Blocking Intensity Index as Wiedenmann et al. (2002)
MGI=apply(totmeridional,c(1,2),mean,na.rm=T)    #Value of meridional gradient inversion
RWB=apply(totrwb,c(1,2),mean,na.rm=T)           #Rossby Wave Breaking parameters
Z500mean=apply(Z500,c(1,2),mean,na.rm=T)          #Z500 mean value

#anticyclonic and cyclonic averages
CN=apply(totrwb,c(1,2),function(x) sum(x[x==(-10)],na.rm=T))/(ndays)*(-10)
ACN=apply(totrwb,c(1,2),function(x) sum(x[x==(10)],na.rm=T))/(ndays)*(10)

print("Saving...")
save(ics,ipsilon,totblocked,frequency,Z500,Z500mean,totmeridional,totBI,totrwb,ACN,CN,MGI,BI,file=outname)


