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
BLOCKDIR=paste(BLOCKDIR,"/",year1,"_",year2,"/",season,"/",sep="")
FIGDIR=paste(FIGDIR,"/Block/",year1,"_",year2,"/",season,"/",sep="")
dir.create(FIGDIR,recursive=T)

#preparing routines
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#which fieds to plot/save
fieldlist=c("InstBlock","Z500","MGI","BI","CN","ACN","BlockEvents","DurationEvents","NumberEvents")

#loading Blocking file
outname=paste(BLOCKDIR,"/Block_",exp,"_",year1,"_",year2,"_",season,sep="")
load(outname)
outname2=paste(BLOCKDIR,"/Events_",exp,"_",year1,"_",year2,"_",season,sep="")
load(outname2)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
namefile=paste(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")


for (var in fieldlist)
{
        #name of the var
        if (var=="InstBlock")
                {longvar="Instantaneous Blocking frequency"; unit="%"; field=frequency}
        if (var=="Z500")
                {longvar="Geopotential Height"; unit="m"; field=Z500mean}
        if (var=="BI")
                {longvar="BI index"; unit=""; field=BI}
        if (var=="MGI")
                {longvar="MGI index"; unit=""; field=MGI}
        if (var=="ACN")
                {longvar="Anticyclonic RWB frequency"; unit="%"; field=ACN}
	if (var=="CN")
                {longvar="Cyclonic RWB frequency"; unit="%"; field=CN}
	if (var=="BlockEvents")
                {longvar="Blocking Events frequency"; unit="%"; field=block$percentage}
	if (var=="DurationEvents")
                {longvar="Blocking Events duration"; unit="days"; field=block$duration}
	if (var=="NumberEvents")
                {longvar="Blocking Events number"; unit=""; field=block$nevents}


        # dimensions definition
        TIME=paste("days since ","1900","-","01","-01 00:00:00",sep="")
        LEVEL=50000
        x <- ncdim_def( "Lon", "degrees", ics)
        y <- ncdim_def( "Lat", "degrees", ipsilon)
        z <- ncdim_def( "Lev", "Pa", LEVEL)
        t <- ncdim_def( "Time", TIME, 1,unlim=T)

        #variable definitions
        var_ncdf=ncvar_def(var,unit,list(x,y,z,t),-999,longname=longvar,prec="single")
        assign(paste("var",var,sep=""),var_ncdf)
        assign(paste("field",var,sep=""),field)
}

namelist=paste("var",fieldlist,sep="")
nclist<-lapply(namelist,get)
ncfile <- nc_create(namefile,nclist)
for (var in fieldlist)
{
        # create ncdf file
        ncvar_put(ncfile, fieldlist[which(var==fieldlist)], get(paste("field",var,sep="")), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
}

nc_close(ncfile)


##########################################################
#-----------------Produce figures------------------------#
##########################################################

#set reference field
dataset_ref="ERAINTERIM"
year1_ref=1979
year2_ref=2014

#open reference field
for (field in fieldlist)
                {
                nomefile=paste(PROGDIR,"/clim/Block/BlockClim_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc",sep="")
                field_ref=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
                assign(paste(field,"_ref",sep=""),field_ref)
                }

#loop on fields
for (field in fieldlist)

	{
	#define field-dependent properties
	if (field=="InstBlock")
	{
		field_exp=frequency; 
		color_field=palette1; color_diff=palette2
		lev_field=seq(0,35,3); lev_diff=seq(-10.5,10.5,1)
		legend_unit="Blocked Days (%)"; title_name="Instantaneous Blocking frequency:"; legend_distance=3
	}

	if (field=="BlockEvents")
	{
		field_exp=block$percentage
                color_field=palette1; color_diff=palette2
                lev_field=seq(0,25,3); lev_diff=seq(-10.5,10.5,1)
                legend_unit="Blocked Days (%)"; title_name="Blocking Events frequency:"; legend_distance=3
	}
	
	if (field=="DurationEvents")
        {
                field_exp=block$duration
                color_field=palette0; color_diff=palette2
                lev_field=seq(5,11.5,.5); lev_diff=seq(-2.1,2.1,.2)
                legend_unit="Duration (days)"; title_name="Duration of Blocking Events:"; legend_distance=3
        }
	
	if (field=="NumberEvents")
        {
                field_exp=block$nevents
                color_field=palette0; color_diff=palette2
                lev_field=seq(0,100,10); lev_diff=seq(-42.5,42.5,5)
                legend_unit=""; title_name="Number of Blocking Events:"; legend_distance=3
        }


	if (field=="Z500")
        {
                field_exp=Z500mean;
                color_field=palette0; color_diff=palette2
                lev_field=seq(4800,6000,50); lev_diff=seq(-310,310,20)
                legend_unit="Geopotential Height (m)"; title_name="Z500:" ; legend_distance=4
        }

	if (field=="BI")
        {
                field_exp=BI;
                color_field=palette0; color_diff=palette2
                lev_field=seq(1,6,0.25); lev_diff=seq(-2.1,2.1,.2)
                legend_unit="BI index"; title_name="Blocking Intensity (BI):" ; legend_distance=3
        }

	if (field=="MGI")
        {
                field_exp=MGI;
                color_field=palette0; color_diff=palette2
                lev_field=seq(0,15,1); lev_diff=seq(-5.25,5.25,.5)
                legend_unit="MGI Index"; title_name="Meridional Gradient Inversion (MGI):" ; legend_distance=3
        }

	if (field=="ACN" | field=="CN")
        {
                if (field=="ACN") {field_exp=ACN; title_name="Anticyclonic Rossby wave breaking frequency:"}
		if (field=="CN") {field_exp=CN; title_name="Cycclonic Rossby wave breaking frequency:"}
                color_field=palette1; color_diff=palette2
                lev_field=seq(0,20,2); lev_diff=seq(-5.25,5.25,.5)
                legend_unit="RWB frequency (%)"; legend_distance=3
        }


	field_ref=get(paste(field,"_ref",sep=""))
	
	
	#secondary plot properties
	nlev_field=length(lev_field)-1
	nlev_diff=length(lev_diff)-1
	lat_lim=c(20,90)
	info_exp=paste(exp,year1,"-",year2,season)
	info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)

	#final plot production
	name=paste(FIGDIR,"/",field,"_",exp,"_",year1,"_",year2,"_",season,".pdf",sep="")
	print(name)
	pdf(file=name,width=12,height=12,onefile=T)
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



