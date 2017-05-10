#!/bin/bash

#program definition

#CDO
#if you CDO is not equipped of NetCDF4 compression change "cdo4" command here
cdo=/usr/bin/cdo
cdonc="$cdo -f nc"
cdo4="$cdo -f nc4 -z zip"

#Rscript is the script-launcher by R
Rscript=/usr/bin/Rscript

#program folder
export PROGDIR=$(pwd)
#data folder
export DATADIR=/work/users/paolo/scratch/miles

####################################
#no need to change below this line #
####################################

#Z500 folder
export ZDIR=$DATADIR/Z500/$exp
#TEMPDIR folder
export TEMPDIR=$DATADIR/tempdir/${exp}_$RANDOM
#NetCDF output dir
export FILESDIR=$DATADIR/files/$exp
#figures folder
export FIGDIR=$DATADIR/figures/$exp

# file type
export output_file_type

#creating folders
mkdir -p $ZDIR $FIGDIR $FILESDIR $TEMPDIR
echo $TEMPDIR

#safety check
echo "Check if R has been loaded"
command -v $Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
Rscript $PROGDIR/config/installpack.R


