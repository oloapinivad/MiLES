#!/bin/bash

#program definition

#CDO
#if you CDO is not equipped of NetCDF4 compression change "cdo4" command here
#Rscript is the script-launcher by R

if [ "$(hostname)" == "lorenzo" ] ; then 
	cdo=/usr/local/apps/cdo/1.8.2/bin/cdo
	Rscript=/usr/local/apps/R/3.1.1/bin/Rscript
else
	cdo=/usr/local/apps/cdo/1.8.2/bin/cdo
	Rscript=/usr/local/apps/R/3.3.1/bin/Rscript
fi

cdonc="$cdo -f nc"
cdo4="$cdo -f nc4 -z zip"

#program folder where MiLES is placed
export PROGDIR=$(pwd)
#data folder where place output (Z500 files, NetCDF files and figures)
export OUTPUTDIR=/scratch/ms/it/ccpd/miles

####################################
#no need to change below this line #
####################################

#NetCDF output dir
export FILESDIR=$OUTPUTDIR/files
#figures folder
export FIGDIR=$OUTPUTDIR/figures

# file type
export output_file_type
export map_projection


#creating folders
mkdir -p $ZDIR $FIGDIR $FILESDIR

#safety check
echo "Check if CDO has been loaded"
command -v cdo -v >/dev/null 2>&1 || { echo "CDO module is not loaded. Aborting." >&2; exit 1; }
echo "CDO found: proceeding..."

echo "Check if R has been loaded"
command -v Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
Rscript $PROGDIR/config/installpack.R


