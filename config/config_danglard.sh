#!/bin/bash

#program definition

#CDO
#program definition
cdo=/opt/local/bin/cdo
cdonc="$cdo -f nc"
#if you CDO is not equipped of NetCDF4 compression change options here
cdo4="$cdo -f nc4 -z zip"

#Rscript is the script-launcher by R
Rscript=/usr/local/bin/Rscript

#program folder where MiLES is placed
export PROGDIR=/Users/paolo/Dropbox/MiLES/MiLES_v0.31
#data folder where place output (Z500 files, NetCDF files and figures)
export OUTPUTDIR=/Users/paolo/Desktop/miles/data

####################################
#no need to change below this line #
####################################

#Z500 folder
export DATADIR=$OUTPUTDIR/Z500
#NetCDF output dir
export FILESDIR=$OUTPUTDIR/files
#figures folder
export FIGDIR=$OUTPUTDIR/figures

# file type
export output_file_type

# if we are using standard climatology
if [[ ${std_clim} -eq 1 ]] ; then
	export dataset_ref="ERAINTERIM"
	export year1_ref=1979
	export year2_ref=2014
	export REFDIR=$PROGDIR/clim
else
	export REFDIR=$FILESDIR
fi


#creating folders
mkdir -p $ZDIR $FIGDIR $FILESDIR $OUTPUTDIR

#safety check
echo "Check if R has been loaded"
command -v $Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
Rscript $PROGDIR/config/installpack.R


