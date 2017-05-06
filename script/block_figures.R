######################################################
#--------Blocking routines plotting for ECE3---------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################


#get environmental variables
PROGDIR<-Sys.getenv(c("PROGDIR"))
BLOCKDIR<-Sys.getenv(c("BLOCKDIR"))
FIGDIR<-Sys.getenv(c("FIGDIRBLOCK"))
#R_LIBLOC<-Sys.getenv(c("R_LIBLOC"))
R_LIBLOC=.libPaths()[1]

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



#loadin packages
library("spam",lib.loc=R_LIBLOC)
library("maps",lib.loc=R_LIBLOC)
library("fields",lib.loc=R_LIBLOC)
library("ncdf4",lib.loc=R_LIBLOC)

#preparing routines
source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#which fieds to plot
fieldlist=c("Block","Z500","BI","MGI")

dataset_ref="ERAINTERIM"; year1_ref=1989; year2_ref=2010


#open reference fields
for (field in fieldlist)
		{
		nomefile=paste(PROGDIR,"/clim/Block/Blocking_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc",sep="")
		field_ref=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
		assign(paste(field,"_ref",sep=""),field_ref)
		}

#color palette to be used
palette0=tim.colors
palette1=colorRampPalette(c("white","orange","darkred"))
palette2=colorRampPalette(c("blue","white","red"))

#loading Blocking file
outname=paste(BLOCKDIR,"/Block_",exp,"_",year1,"_",year2,"_",season,sep="")
load(outname)

#loop on fields
for (field in fieldlist)

	{
	#define field-dependent properties
	if (field=="Block")
	{
		field_exp=frequency; 
		color_field=palette1; color_diff=palette2
		lev_field=seq(0,35,3); lev_diff=seq(-10.5,10.5,1)
		legend_unit="Blocked Days (%)"; title_name="Instantaneous Blocking:"; legend_distance=3
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
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_field),max(lev_field)),col=color_field(nlev_field),legend.args=list(side=4,line=legend_distance,cex=1.2,text=legend_unit))

	#delta field plot
	filled.contour3(ics,ipsilon,field_exp-field_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=color_diff,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_diff),max(lev_diff)),col=color_diff(nlev_diff),legend.args=list(side=4,line=legend_distance,cex=1.2,text=legend_unit))

	dev.off()
	}

#saving output to netcdf files
print("saving NetCDF climatologies...")
filein=paste(BLOCKDIR,"/Block_",exp,"_",year1,"_",year2,"_",season,sep="")
load(filein)
namefile=paste(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")

# create ncdf file and define variables to put in
varlist=c("Block","Z500","MGI","BI")

for (var in varlist)
{
        #name of the var
        if (var=="Block")
                {longvar="Instantaneous Blocking"; unit="% of Blocked Days"; field=frequency}
        if (var=="Z500")
                {longvar="Geopotential Height"; unit="m"; field=Z500mean}
        if (var=="BI")
                {longvar="BI index"; unit=""; field=BI}
         if (var=="MGI")
                {longvar="MGI index"; unit=""; field=MGI}

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

namelist=paste("var",varlist,sep="")
nclist<-mget(namelist)
ncfile <- nc_create(namefile,nclist)
for (var in varlist)
{
        # create ncdf file
        ncvar_put(ncfile, varlist[which(var==varlist)], get(paste("field",var,sep="")), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
}

nc_close(ncfile)


