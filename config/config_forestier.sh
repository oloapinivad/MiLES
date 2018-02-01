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
export PROGDIR=/Users/paolo/Desktop/MiLES
#data folder where place output (Z500 files, NetCDF files and figures)
export OUTPUTDIR=/Users/paolo/Desktop/milesdata

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

