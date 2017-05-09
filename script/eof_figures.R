######################################################
#--------Routines for EOFs plotting for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################


#get environmental variables
PROGDIR<-Sys.getenv(c("PROGDIR"))
EOFDIR0<-Sys.getenv(c("EOFDIR"))
FIGDIR0<-Sys.getenv(c("FIGDIREOFS"))

#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]
tele=args[5]

#correct folder to experiment dependent
FIGDIR=paste(FIGDIR0,"/EOFs/",tele,"/",year1,"_",year2,"/",season,"/",sep="")
EOFDIR=paste(EOFDIR0,"/",tele,"/",year1,"_",year2,"/",season,"/",sep="")
dir.create(FIGDIR,recursive=T)

#preparing routines
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#EOFs to plot (depends on how many computed by CDO!)
neofs=4

#loading anomalies and variances of experiment
print(EOFDIR)
nomefile=paste(EOFDIR,"Z500_monthly_anomalies_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
anomalies_exp=ncdf.opener(nomefile,"zg","lon","lat",rotate=T)
nomefile=paste(EOFDIR,tele,"_Z500_eigenvalues_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
variance=ncdf.opener(nomefile,"zg")
variance_exp=round(variance[neofs]/sum(variance)*100,1)


#set reference field
dataset_ref="ERAINTERIM"
year1_ref=1989
year2_ref=2010

#loading anomalies and variance of reference
info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)

if (dataset_ref=="ERAINTERIM" & year1_ref=="1989" & year2_ref=="2010")
	{refdir=paste(PROGDIR,"/clim/EOFs/",tele,"/",season,"/",sep="")} else {refdir=paste0(gsub(exp,dataset_ref,EOFDIR0),"/",tele,"/",year1_ref,"_",year2_ref,"/",season,"/")}

nomefile=paste(refdir,"Z500_monthly_anomalies_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc",sep="")
anomalies_ref=ncdf.opener(nomefile,"zg","lon","lat",rotate=T)
nomefile=paste(refdir,tele,"_Z500_eigenvalues_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc",sep="")
variance=ncdf.opener(nomefile,"zg")
variance_ref=round(variance[1:neofs]/sum(variance)*100,1)

#loop on number of EOFs
for (neof in 1:neofs)
	{

	#loading PCs of experiment and normalize
	nomefile=paste(EOFDIR,"/",tele,"_monthly_timeseries_",exp,"_",year1,"_",year2,"_",season,"_0000",neof-1,".nc",sep="")
	timeseries_exp=ncdf.opener(nomefile,"zg")
	timeseries_exp=(timeseries_exp-mean(timeseries_exp))/sd(timeseries_exp)


	#linear regression on Z500 anomalies for experiment (faster function)
	#linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lm(linreg ~ timeseries_exp,na.action=na.exclude)$coef[2])
	linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_exp,ncol=1),linreg)$coefficients)

	#loading PC of reference and normalize
	nomefile=paste(refdir,tele,"_monthly_timeseries_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,"_0000",neof-1,".nc",sep="")
	timeseries_ref=ncdf.opener(nomefile,"zg")
	timeseries_ref=(timeseries_ref-mean(timeseries_ref))/sd(timeseries_ref)

	#linear regression on Z500 anomalies for reference (faster lm.fit function)
	#linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lm(linreg ~ timeseries_ref,na.action=na.exclude)$coef[2])
	linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_ref,ncol=1),linreg)$coefficients)

	#check and flip signs
	if (cor(c(linear_ref),c(linear_exp))<0) {linear_exp=(-linear_exp)}

	#-----plotting-------#
	
	#plot properties
	lev_field=seq(-150,150,20)
	lev_diff=seq(-95,95,10)
	nlev_field=length(lev_field)-1
	nlev_diff=length(lev_diff)-1
	lat_lim=c(20,90)
	title_name=paste("EOF",neof,sep="")
	info_exp=paste(exp,year1,"-",year2,season)

	#final plot production
	figname=paste(FIGDIR,"EOF",neof,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type,sep="")
	print(figname)

	# Chose output format for figure - by JvH
        if (tolower(output_file_type) == "png") {
           png(filename = figname, width=png_width, height=png_height)
        } else if (tolower(output_file_type) == "pdf") {
            pdf(file=figname,width=pdf_width,height=pdf_height,onefile=T)
        } else if (tolower(output_file_type) == "eps") {
            setEPS(width=pdf_width,height=pdf_height,onefile=T,paper="special")
            postscript(figname)
        }

	#plot properties
	par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))
	print(quantile(linear_exp))

	filled.contour3(ics,ipsilon,linear_exp,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_exp),levels=lev_field,color.palette=palette0,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	text(120,85,paste("Variance Explained: ",variance_exp[neof],"%",sep=""),cex=2)

	filled.contour3(ics,ipsilon,linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_ref),levels=lev_field,color.palette=palette0,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_field,color.palette=palette0,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)
	text(120,85 ,paste("Variance Explained: ",variance_ref[neof],"%",sep=""),cex=2)

	#delta field plot
	filled.contour3(ics,ipsilon,linear_exp-linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=palette2,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_field,color.palette=palette2,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)
	
	dev.off()
	}
