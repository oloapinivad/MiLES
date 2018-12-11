rm(list=ls())
PROGDIR="/home/paolo/MiLES"
FILESDIR="/work/users/paolo/miles/files"
source(file.path(PROGDIR,"script/basis_functions.R"))
dataset="ERAI"
project=expid=ens=NA
year1=1979
year2=1980
season="DJF"
z500filename="/work/users/paolo/miles/data/zg500/ERAI/zg500_ERAI_fullfile.nc"
#z500filename="/work/users/paolo/miles/data/zg500/CMIP5/bcc-csm1-1/historical/r1/zg500_bcc-csm1-1_historical_r1_fullfile.nc"

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#load file
fieldlist=ncdf.opener.universal(z500filename,namevar="zg",tmonths=timeseason,tyears=years,rotate="full")
print(str(fieldlist))

#extract calendar and time unit from the original file
timeaxis=fieldlist$time

#time array to simplify time filtering
etime=power.date.new(timeaxis,verbose=T)


#declare variable
Z500=fieldlist$field
Z500mean=apply(Z500,c(1,2),mean)

# list of isohypses on which evaluate the MI
isolvls=seq(4900,6200,5)

# reference latitude (60N following Di Capua et al., 2016)
ref_lat=60

savefile1=file.builder(FILESDIR,"Block","BlockClim",project,dataset,expid,ens,year1,year2,season)


#saving output to netcdf files
print("saving NetCDF climatologies...")

#which fieds to plot/save
full_savelist=savelist="Z500"

# dimension definition (using default 1850-01-01 reftime)
dims=ncdf.defdims(ics,ipsilon,timeaxis)

nc_field=nc_var<-vector("list",length(savelist))
nc_fullfield=nc_fullvar<-vector("list",length(full_savelist))

for (var in savelist)
{
        #name of the var
        if (var=="Z500")
                {longvar="Geopotential Height"; unit="m"; field=Z500mean; full_field=Z500}

                var_ncdf=ncvar_def(var,unit,list(dims$x,dims$y,dims$z,t=dims$tclim),-999,longname=longvar,prec="single",compression=1)
                full_var_ncdf=ncvar_def(var,unit,list(dims$x,dims$y,dims$z,t=dims$t),-999,longname=longvar,prec="single",compression=1)
	 #assign(paste0("var",var),var_ncdf)
        #assign(paste0("full_var",var),full_var_ncdf)
        nc_var[[which(var==savelist)]] <- var_ncdf    
        nc_field[[which(var==savelist)]] <- field

        #assign(paste0("field",var),field)
        #assign(paste0("full_field",var),full_field)
        nc_fullvar[[which(var==savelist)]] <- var_ncdf
        nc_fullfield[[which(var==savelist)]] <- field

}

#ncdf.write(savefile2,full_savelist,vartype="full")

