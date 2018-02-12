#!/bin/bash

####################################
# ---- folder definition --------- #
####################################

# INDIR ->data folder: all the geopotential height data should be here
# you must customozie this according to the dataset you analyze and your local file structure
# dataset and ens is called by the z500_prepare script
if [ "${dataset}" == "ERA5" ] |  [ "${dataset}" == "ERAI" ] ; then
    INDIR=/scratch/rd/nedd/regimes/${dataset}_daily
    expected_input_name=*${dataset}*${ens}.grb
fi

if [ "${dataset}" == "S4" ] ; then
    INDIR=/home/ms/it/ccpd/regimes/S4_daily_nov
    expected_input_name=*/*S4*_${ens}.grb
fi

if [ "${dataset}" == "S5" ] |[ "${dataset}" == "S3" ] ; then
    INDIR=/scratch/rd/nedd/regimes/${dataset}_24h_nov
    expected_input_name=*/*${dataset}*_${ens}.grb
fi

if [ "${dataset}" == "S4AMIP" ] | [ "${dataset}" == "S5AMIP" ] | [ "${dataset}" == "S5AMIP_ERAI" ] ; then
    INDIR=/gpfs/scratch/ms/it/ccpd/mars/${dataset}
    expected_input_name=*/*S0*_${ens}.grb
fi


if [ "${dataset}" == "S5LR" ] ; then
    INDIR=/scratch/rd/nedd/regimes/S5_24h_nov_lowres
    expected_input_name=*/*S5*_${ens}.grb
fi

#program folder where MiLES is placed
export PROGDIR=$(pwd)
#data folder where place output (Z500 files, NetCDF files and figures)
export OUTPUTDIR=/scratch/ms/it/ccpd/miles

####################################
# ----  program definition  ------ #
####################################

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
command -v $cdo -v >/dev/null 2>&1 || { echo "CDO module is not loaded. Aborting." >&2; exit 1; }
echo "CDO found: proceeding..."

echo "Check if R has been loaded"
command -v $Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

#echo "Check if NetCDF  has been loaded"
#command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
#echo "NetCDF found: starting...."

#R check for key packages
$Rscript $PROGDIR/config/installpack.R


