######################################################
#--------Routines for EOFs plotting for ECE3---------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################


#get environmental variables
PROGDIR<-Sys.getenv(c("PROGDIR"))
EOFDIR<-Sys.getenv(c("EOFDIR"))
FIGDIR<-Sys.getenv(c("FIGDIREOFS"))
#R_LIBLOC<-Sys.getenv(c("R_LIBLOC"))
R_LIBLOC=.libPaths()[1]

#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]
tele=args[5]

#correct folder to experiment dependent
FIGDIR=paste(FIGDIR,"/EOFs/",tele,"/",year1,"_",year2,"/",season,"/",sep="")
EOFDIR=paste(EOFDIR,"/",tele,"/",year1,"_",year2,"/",season,"/",sep="")
dir.create(FIGDIR,recursive=T)

#loadin packages
library("spam",lib.loc=R_LIBLOC)
library("maps",lib.loc=R_LIBLOC)
library("fields",lib.loc=R_LIBLOC)
library("ncdf4",lib.loc=R_LIBLOC)

#preparing routines
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#loading anomalies and variances of experiment
print(EOFDIR)
nomefile=paste(EOFDIR,"Z500_monthly_anomalies_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
anomalies_exp=ncdf.opener(nomefile,"zg","lon","lat",rotate=T)
nomefile=paste(EOFDIR,tele,"_Z500_eigenvalues_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
variance=ncdf.opener(nomefile,"zg")
variance_exp=round(variance[1:4]/sum(variance)*100,1)

#loading anomalies and variance of reference
info_ref=paste("ERAINTERIM 1989 - 2010",season)
refdir=paste(PROGDIR,"/clim/EOFs/",tele,"/",season,"/",sep="")
nomefile=paste(refdir,"Z500_monthly_anomalies_ERAINTERIM_1989_2010_",season,".nc",sep="")
anomalies_ref=ncdf.opener(nomefile,"zg","lon","lat",rotate=T)
nomefile=paste(refdir,tele,"_Z500_eigenvalues_ERAINTERIM_1989_2010_",season,".nc",sep="")
variance=ncdf.opener(nomefile,"zg")
variance_ref=round(variance[1:4]/sum(variance)*100,1)

#loop on number of EOFs
for (neof in c(1,2,3,4))
	{

	#loading PCs of experiment
	nomefile=paste(EOFDIR,"/",tele,"_monthly_timeseries_",exp,"_",year1,"_",year2,"_",season,"_0000",neof-1,".nc",sep="")
	timeseries_exp=ncdf.opener(nomefile,"zg")
	timeseries_exp=(timeseries_exp-mean(timeseries_exp))/sd(timeseries_exp)


	#linear regression on Z500 anomalies for experiment
	linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lm(linreg ~ timeseries_exp,na.action=na.exclude)$coef[2])

	nomefile=paste(refdir,tele,"_monthly_timeseries_ERAINTERIM_1989_2010_",season,"_0000",neof-1,".nc",sep="")
	timeseries_ref=ncdf.opener(nomefile,"zg")
	timeseries_ref=(timeseries_ref-mean(timeseries_ref))/sd(timeseries_ref)

	#linear regression on Z500 anomalies for reference
	linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lm(linreg ~ timeseries_ref,na.action=na.exclude)$coef[2])

	#check and flip signs
	if (cor(c(linear_ref),c(linear_exp))<0) {linear_exp=(-linear_exp)}

	#-----plotting-------#
	
	#color palette to be used
	palette0=tim.colors
	palette1=colorRampPalette(c("white","orange","darkred"))
	palette2=colorRampPalette(c("blue","white","red"))

	#plot properties
	lev_field=seq(-150,150,20)
	lev_diff=seq(-95,95,10)
	nlev_field=length(lev_field)-1
	nlev_diff=length(lev_diff)-1
	lat_lim=c(20,90)
	title_name=paste("EOF",neof,sep="")
	info_exp=paste(exp,year1,"-",year2,season)

	#final plot production
	name=paste(FIGDIR,"EOF",neof,"_",exp,"_",year1,"_",year2,"_",season,".pdf",sep="")
	print(name)
	pdf(file=name,width=12,height=12,onefile=T)
	par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))
	print(quantile(linear_exp))

	filled.contour3(ics,ipsilon,linear_exp,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_exp),levels=lev_field,color.palette=palette0,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	text(140,80,paste("Variance Explained: ",variance_exp[neof],"%",sep=""),cex=1.4)

	filled.contour3(ics,ipsilon,linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_ref),levels=lev_field,color.palette=palette0,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_field),max(lev_field)),col=palette0(nlev_field),legend.args=list(side=4,line=3,cex=1.2,text="m"))
	text(140,80,paste("Variance Explained: ",variance_ref[neof],"%",sep=""),cex=1.4)

	#delta field plot
	filled.contour3(ics,ipsilon,linear_exp-linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=palette2,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_diff),max(lev_diff)),col=palette2(nlev_diff),legend.args=list(side=4,line=3,cex=1.2,text="m"))

	
	dev.off()
	}
