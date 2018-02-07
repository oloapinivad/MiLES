######################################################
#--------Routines for EOFs plotting for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

miles.eof.figures<-function(exp,ens,year1,year2,dataset_ref,ens_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,tele)
{

#R configuration file 
source(CFGSCRIPT)

#use filebuilding script to access to file
nomefile_exp=file.builder(FILESDIR,paste0("EOFs_beta/",tele),"EOFs",exp,ens,year1,year2,season)

# check for REFDIR==FILESDIR, i.e. if we are using the climatology provided by MiLES or another dataset MiLES-generated 
    if (REFDIR!=FILESDIR) {
        nomefile_ref=paste0(file.path(REFDIR,paste0("EOFs/",tele)),"/EOFs_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
    } else {

        #use file.builder to create the path of the blocking files
        nomefile_ref=file.builder(FILESDIR,paste0("EOFs_beta/",tele),"/EOFs_",dataset_ref,ens_ref,year1_ref,year2_ref,season)
    }

#EOFs to plot (depends on how many computed by CDO!)
neofs=4

##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

#loading anomalies and variances of experiment
anomalies_exp=ncdf.opener(nomefile_anomalies,"zg","lon","lat",rotate="full")
variance=ncdf.opener(nomefile_variance,"zg")
variance_exp=round(variance[1:neofs]/sum(variance)*100,1)

#loading reference field
anomalies_ref=ncdf.opener(nomefile_ref_anomalies,"zg","lon","lat",rotate="full")
variance=ncdf.opener(nomefile_ref_variance,"zg")
variance_ref=round(variance[1:neofs]/sum(variance)*100,1)


##########################################################
#-----------------Produce figures------------------------#
##########################################################

#plot properties
if (ens=="NO") {info_exp=paste(exp,year1,"-",year2,season)} else {info_exp=paste(exp,ens,year1,"-",year2,season)}
if (ens_ref=="NO") {info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)} else {info_ref=paste(dataset_ref,ens_ref,year1_ref,"-",year2_ref,season)}
lev_field=seq(-150,150,20)
lev_diff=seq(-95,95,10)

#loop on number of EOFs
for (neof in 1:neofs) {

	#loading PCs of experiment and normalize
    if (ens=="NO") {
        nomefile=paste(EOFDIR,"/",tele,"_monthly_timeseries_",exp,"_",year1,"_",year2,"_",season,"_0000",neof-1,".nc",sep="")
    } else {
        nomefile=paste(EOFDIR,"/",tele,"_monthly_timeseries_",exp,"_",ens,"_",year1,"_",year2,"_",season,"_0000",neof-1,".nc",sep="")
    }
	timeseries_exp0=ncdf.opener(nomefile,"zg")
	timeseries_exp=standardize(timeseries_exp0)

	#linear regression on Z500 anomalies for experiment (faster function)
	#linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lm(linreg ~ timeseries_exp,na.action=na.exclude)$coef[2])
	linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_exp,ncol=1),linreg)$coefficients)
	
	#loading PC of reference and normalize
    if (ens=="NO") {
	    nomefile=paste(REFDIR,"/",tele,"_monthly_timeseries_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,"_0000",neof-1,".nc",sep="")
    } else {
        nomefile=paste(REFDIR,"/",tele,"_monthly_timeseries_",dataset_ref,"_",ens_ref,"_",year1_ref,"_",year2_ref,"_",season,"_0000",neof-1,".nc",sep="")
    }
	timeseries_ref0=ncdf.opener(nomefile,"zg")
	timeseries_ref=standardize(timeseries_ref0)

	#linear regression on Z500 anomalies for reference (faster lm.fit function)
	#linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lm(linreg ~ timeseries_ref,na.action=na.exclude)$coef[2])
	linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_ref,ncol=1),linreg)$coefficients)

	#check and flip signs (to be in agreement with reference field)
	if (cor(c(linear_ref),c(linear_exp))<0) {linear_exp=(-linear_exp)}
	
	#-----plotting-------#
	
	#plot properties
	region=tele
	if (tele=="NAO") {region="North Atlantic"}
	if (tele=="AO") {region="Northern Hemisphere"}
	title_name=paste0(region," EOF",neof)

	#final plot production
    if (ens=="NO") {
	    figname=paste0(FIGDIREOF,"/EOF",neof,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type)
    } else {
        figname=paste0(FIGDIREOF,"/EOF",neof,"_",exp,"_",ens,"_",year1,"_",year2,"_",season,".",output_file_type)
    }
	print(figname)

	# Chose output format for figure - by JvH
    open.plot.device(figname,output_file_type,CFGSCRIPT)

	#plot properties
	par(plotpar)

	im=plot.prepare(ics,ipsilon,linear_exp,proj=map_projection,lat_lim=lat_lim)
        filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_exp),levels=lev_field,color.palette=palette3,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
        mtext(title_name,side=3,line=.5,outer=TRUE,cex=2,font=2)
        proj.addland(proj=map_projection)
	text(120,85,paste("Variance Explained: ",variance_exp[neof],"%",sep=""),cex=2)

	im=plot.prepare(ics,ipsilon,linear_ref,proj=map_projection,lat_lim=lat_lim)
        filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_exp),levels=lev_field,color.palette=palette3,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
        mtext(title_name,side=3,line=.5,outer=TRUE,cex=2,font=2)
        proj.addland(proj=map_projection)
	image.scale3(volcano,levels=lev_field,color.palette=palette3,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1*af,line.label=3)
	text(120,85 ,paste("Variance Explained: ",variance_ref[neof],"%",sep=""),cex=2)

	#delta field plot
        im=plot.prepare(ics,ipsilon,linear_exp-linear_ref,proj=map_projection,lat_lim=lat_lim)        
	filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste("Difference"),levels=lev_diff,color.palette=palette2,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
        proj.addland(proj=map_projection)
        image.scale3(volcano,levels=lev_diff,color.palette=palette2,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1*af,line.label=3)

	
	dev.off()
	}

}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","ens","year1","year2","dataset_ref","ens_ref","year1_ref","year2_ref","season","FIGDIR","FILESDIR","REFDIR","CFGSCRIPT","PROGDIR","tele")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
        print(paste("If running from bash, please specify the",req_args,"arguments here below:"))
        print(name_args)
     } else {
# when the number of arguments is ok run the function()
        for (k in 1:req_args) {assign(name_args[k],args[k])}
        source(paste0(PROGDIR,"/script/basis_functions.R"))
        miles.eof.figures(exp,ens,year1,year2,dataset_ref,ens_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,tele)
     }
}


