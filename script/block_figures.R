######################################################
#--Blocking routines plotting and saving for MiLES---#
#-------------P. Davini (Oct 2014)-------------------#
######################################################

#get environmental variables
PROGDIR<-Sys.getenv(c("PROGDIR"))
BLOCKDIR<-Sys.getenv(c("BLOCKDIR"))
FIGDIR<-Sys.getenv(c("FIGDIRBLOCK"))

#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]

#correct folder to year and season dependence
BLOCKDIR=paste(BLOCKDIR,"/",year1,"_",year2,"/",season,sep="")
FIGDIR=paste(FIGDIR,"/Block/",year1,"_",year2,"/",season,sep="")
dir.create(FIGDIR,recursive=T)

#preparing routines
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#which fieds to plot/save
fieldlist=c("InstBlock","Z500","MGI","BI","CN","ACN","BlockEvents","DurationEvents","NumberEvents")

##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

#open reference field
for (field in fieldlist)
                {
                nomefile=paste(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
                field_exp=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
                assign(paste(field,"_exp",sep=""),field_exp)
                }

#set reference field
dataset_ref="ERAINTERIM"
year1_ref=1979
year2_ref=2014

#open reference field
for (field in fieldlist)
                {
		if (dataset_ref=="ERAINTERIM" & year1_ref=="1979" & year2_ref=="2014")
			{
			nomefile=paste0(PROGDIR,"/clim/Block/BlockClim_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")} else {
			nomefile=paste0(gsub(exp,dataset_ref,BLOCKDIR),"/BlockClim_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
			}
                field_ref=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
                assign(paste(field,"_ref",sep=""),field_ref)
                }

##########################################################
#-----------------Produce figures------------------------#
##########################################################

#loop on fields
for (field in fieldlist)

	{
	#define field-dependent properties
	if (field=="InstBlock")
	{
		color_field=palette1; color_diff=palette2
		lev_field=seq(0,35,3); lev_diff=seq(-10.5,10.5,1)
		legend_unit="Blocked Days (%)"; title_name="Instantaneous Blocking frequency:"; legend_distance=3
	}

	if (field=="BlockEvents")
	{
                color_field=palette1; color_diff=palette2
                lev_field=seq(0,25,3); lev_diff=seq(-10.5,10.5,1)
                legend_unit="Blocked Days (%)"; title_name="Blocking Events frequency:"; legend_distance=3
	}
	
	if (field=="DurationEvents")
        {
                color_field=palette0; color_diff=palette2
                lev_field=seq(5,11.5,.5); lev_diff=seq(-2.1,2.1,.2)
                legend_unit="Duration (days)"; title_name="Duration of Blocking Events:"; legend_distance=3
        }
	
	if (field=="NumberEvents")
        {
                color_field=palette0; color_diff=palette2
                lev_field=seq(0,100,10); lev_diff=seq(-42.5,42.5,5)
                legend_unit=""; title_name="Number of Blocking Events:"; legend_distance=3
        }

	if (field=="Z500")
        {
                color_field=palette0; color_diff=palette2
                lev_field=seq(4800,6000,50); lev_diff=seq(-310,310,20)
                legend_unit="Geopotential Height (m)"; title_name="Z500:" ; legend_distance=4
        }

	if (field=="BI")
        {
                color_field=palette0; color_diff=palette2
                lev_field=seq(1,6,0.25); lev_diff=seq(-2.1,2.1,.2)
                legend_unit="BI index"; title_name="Blocking Intensity (BI):" ; legend_distance=3
        }

	if (field=="MGI")
        {
                color_field=palette0; color_diff=palette2
                lev_field=seq(0,15,1); lev_diff=seq(-5.25,5.25,.5)
                legend_unit="MGI Index"; title_name="Meridional Gradient Inversion (MGI):" ; legend_distance=3
        }

	if (field=="ACN" | field=="CN")
        {
                if (field=="ACN") {title_name="Anticyclonic Rossby wave breaking frequency:"}
		if (field=="CN") {title_name="Cyclonic Rossby wave breaking frequency:"}
                color_field=palette1; color_diff=palette2
                lev_field=seq(0,20,2); lev_diff=seq(-5.25,5.25,.5)
                legend_unit="RWB frequency (%)"; legend_distance=3
        }

	field_ref=get(paste(field,"_ref",sep=""))
	field_exp=get(paste(field,"_exp",sep=""))
	
	
	#secondary plot properties
	nlev_field=length(lev_field)-1
	nlev_diff=length(lev_diff)-1
	lat_lim=c(20,90)
	info_exp=paste(exp,year1,"-",year2,season)
	info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)

	#final plot production
	figname=paste(FIGDIR,"/",field,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type,sep="")
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

	#panels option
	par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))

	#main experiment plot
	filled.contour3(ics,ipsilon,field_exp,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_exp),levels=lev_field,color.palette=color_field,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)

	#reference field plot
	filled.contour3(ics,ipsilon,field_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_ref),levels=lev_field,color.palette=color_field,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_field,color.palette=color_field,colorbar.label=legend_unit,cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=legend_distance,line.colorbar=1.5)

	#delta field plot
	filled.contour3(ics,ipsilon,field_exp-field_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=color_diff,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_diff,color.palette=color_diff,colorbar.label=legend_unit,cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=legend_distance,line.colorbar=1.5)

	dev.off()
	}



