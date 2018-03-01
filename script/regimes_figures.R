######################################################
#------Regimes routines figures for MiLES------------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

miles.regimes.figures<-function(exp,ens,year1,year2,dataset_ref,ens_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,nclusters)
{

if (nclusters!=4 | season!="DJF") {stop("Beta version: unsupported season and/or number of clusters")}

#R configuration file 
source(CFGSCRIPT)


##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

# loading anomalies and variances of experiment
nomefile=file.builder(FILESDIR,"Regimes","RegimesPattern",exp,ens,year1,year2,season)
frequencies_exp=ncdf.opener(nomefile,"Frequencies")
regimes_exp=ncdf.opener(nomefile,"Regimes","Lon","Lat",rotate="no")

# loading reference field
# check for REFDIR==FILESDIR, i.e. if we are using the climatology provided by MiLES or another dataset MiLES-generated 
    if (REFDIR!=FILESDIR) {
        nomefile_ref=paste0(file.path(REFDIR,"Regimes"),"/RegimesPattern_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
    } else {

        #use file.builder to create the path of the blocking files
        nomefile_ref=file.builder(FILESDIR,"Regimes","RegimesPattern",dataset_ref,ens_ref,year1_ref,year2_ref,season)
    }

#nomefile=paste0(REFDIR,"/RegimesPattern_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
frequencies_ref=ncdf.opener(nomefile_ref,"Frequencies")
regimes_ref=ncdf.opener(nomefile_ref,"Regimes","Lon","Lat",rotate="no")

#try to assign the 4 standard regimes names to the dataset using the distance between 
#the location of the maximum/minimum of the pattern and 4 "standard" locations
#when something is wrong (i.e. multiple assignments) general "Regime X" names are set
#It is not perfect, it is just aimed at simplying the plots 
for (ii in c("ref","exp"))
{
	compose=get(paste0("regimes_",ii))
	names=paste("Regimes",1:nclusters)
    position=rbind(c(-50,60),c(-40,50),c(0,60),c(-15,60))
    rownames(position)<-c("NAO-","Atlantic Ridge","Scandinavian Blocking","NAO+")
    for (i in 1:nclusters)
                {
                MM=which(compose[,,i]==max(compose[,,i],na.rm=T),arr.ind=T)
                mm=which(compose[,,i]==min(compose[,,i],na.rm=T),arr.ind=T)
                if (max(compose[,,i],na.rm=T)>abs(min(compose[,,i],na.rm=T)))
                        {
                        distMM=dist(rbind(c(ics[MM[1]],ipsilon[MM[2]]),position))
                        } else {
                        distMM=dist(rbind(c(ics[mm[1]],ipsilon[mm[2]]),position))
                        }
                #print(distMM)
                names[i]=rownames(position)[which.min(distMM[1:nclusters])]
		if (i>1 & any(names[i]==names[1:max(c(1,i-1))])) {names[i]=paste("Regime",i)	}
		#print(names[i])
                }

assign(paste0("names_",ii),names)
}

#plot properties
lev_field=seq(-250,250,20)
lev_diff=seq(-150,150,20)
#standard properties
if (ens=="NO") {info_exp=paste(exp,year1,"-",year2,season)} else {info_exp=paste(exp,ens,year1,"-",year2,season)}
if (ens_ref=="NO") {info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)} else {info_ref=paste(dataset_ref,ens_ref,year1_ref,"-",year2_ref,season)}

kk0=1
# loop on regimes
for (name in names_ref)
{
	#-----plotting-------#

	#a bit complicated but it is used to compare similar regimes even if they not
	#equal percentage of occurrence.
	ii=which(name==names_exp)
	jj=which(name==names_ref)
	if (length(ii)==0) {ii=which(setdiff(names_exp,names_ref)[kk0]==names_exp); kk0=kk0+1}
	print(name)

    #final plot production
    figname=fig.builder(FIGDIR,"Regimes",paste0("Regime",ii),exp,ens,year1,year2,season,output_file_type)
    print(figname)
    
    # Chose output format for figure - by JvH
    open.plot.device(figname,output_file_type,CFGSCRIPT)
    
    # where to plot frequencies values
    if (map_projection=="no") {varpoints=c(120,85)} else {varpoints=c(0,0.7)}

    #plot properties
	par(plotpar)

	im=plot.prepare(ics,ipsilon,regimes_exp[,,ii],proj=map_projection,lat_lim=lat_lim)
	filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_exp),levels=lev_field,color.palette=palette3,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
        mtext(name,side=3,line=.5,outer=TRUE,cex=2,font=2)
        proj.addland(proj=map_projection)
        text(varpoints[1],varpoints[2],paste("Frequencies: ",round(frequencies_exp[ii],2),"%",sep=""),cex=2)

	im=plot.prepare(ics,ipsilon,regimes_ref[,,ii],proj=map_projection,lat_lim=lat_lim)
        filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_ref),levels=lev_field,color.palette=palette3,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
        proj.addland(proj=map_projection)
        text(varpoints[1],varpoints[2],paste("Frequencies: ",round(frequencies_ref[ii],2),"%",sep=""),cex=2)
        image.scale3(volcano,levels=lev_field,color.palette=palette3,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1*af,line.label=3)

        #delta field plot
	im=plot.prepare(ics,ipsilon,regimes_exp[,,ii]-regimes_ref[,,jj],proj=map_projection,lat_lim=lat_lim)
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
name_args=c("exp","ens","year1","year2","dataset_ref","ens_ref","year1_ref","year2_ref","season","FIGDIR","FILESDIR","REFDIR","CFGSCRIPT","PROGDIR","nclusters")
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
        miles.regimes.figures(exp,ens,year1,year2,dataset_ref,ens_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,nclusters)
    }
}



