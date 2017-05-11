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
#NetCDF output dir
export FILESDIR=$DATADIR/files
#figures folder
export FIGDIR=$DATADIR/figures

# file type
export output_file_type

# if we are using standard climatology
if [[ ${std_clim} -eq 1 ]] ; then
	export dataset_ref="ERAINTERIM"
	export year1_ref=1979
	export year2_ref=2014
	export REFDIR=$PROGDIR/clim/Block/
else
	export REFDIR=$FILESDIR
fi


#creating folders
mkdir -p $ZDIR $FIGDIR $FILESDIR

#safety check
echo "Check if R has been loaded"
command -v $Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
Rscript $PROGDIR/config/installpack.R


