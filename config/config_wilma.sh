#!/bin/bash

####################################
# ---- folder definition --------- #
####################################


# INDIR ->data folder: all the input dataa should be here
# you must customozie this according to the dataset you analyze and your local file structure
# you can use $dataset, $expid and $ens variables for data structure
# also $varname and $level can be used
# the $project variable is an extra element for directory tree
# The combination of the three variables is the standard for CMIP-like experiments.
if [[ $varname == "zg" ]] ; then findvar=Z500 ; fi
if [[ $varname == "ua" ]] ; then findvar=U500 ; fi

if [[ ${project} == "CMIP5" ]] ; then INDIR=$WORK/data/${project}/${dataset}/${expid}/${ens}/day/${findvar} ; fi
if [[ ${project} == "CMIP3" ]] ; then INDIR=$ARCHIVE/work/${project}/${dataset}/${expid}/${ens}/day/${findvar} ; fi 
if [[ "${dataset}" == "NCEP" ]] && [[ $varname == "zg" ]]  ;  then INDIR=$WORK/data/${dataset}/day/hgt ; fi
if [[ "${dataset}" == "ERA40" ]] || [[ "${dataset}" == "ERAI"  ]] ; then INDIR=$WORK/data/${dataset}/day/${findvar} ; fi
if [[ "${dataset}" == "ERAI"  ]] && [[ $varname == "ua" ]] ; then INDIR=/work/datasets/obs/ERAINT/v0/data/${findvar}/6hrs ; fi
if [[ "${dataset}" == "20CRv2c" ]] ;  then INDIR=$WORK/data/${dataset}/${ens}/day/${findvar} ; fi
if [[ "${dataset}" == "CMCC-CM2" ]] ;  then INDIR=$SCRATCH/cmcc ; fi
if [[ "${project}" == "SPHINX" ]] ;  then INDIR=$WORK/data/${project}/${ens}/day/${findvar} ; fi

# to look for some specific file structure
# if commented the program will look for all the netcdf or grib files in the folder
#expected_input_name=*.nc

#program folder where MiLES is placed
export PROGDIR=$HOME/MiLES
#data folder where place output (Z500 files, NetCDF files and figures)
export OUTPUTDIR=$WORK/miles

####################################
# ----  program definition  ------ #
####################################

#program definition

#CDO
#if you CDO is not equipped of NetCDF4 compression change "cdo4" command here
cdo=/usr/bin/cdo
cdonc="$cdo -f nc"
cdo4="$cdo -f nc4 -z zip"

#Rscript is the script-launcher by R
Rscript=/usr/bin/Rscript


####################################
#no need to change below this line #
####################################

#NetCDF output dir
export FILESDIR=$OUTPUTDIR/files
#figures folder
export FIGDIR=$OUTPUTDIR/figures

# file type and export default
export output_file_type=${output_file_type:-pdf}
export map_projection=${map_projection:-azequalarea}
export expected_file_output=${expected_file_output:-}

#creating folders
mkdir -p $FIGDIR $FILESDIR

#safety check
echo "Check if R has been loaded"
command -v $Rscript >/dev/null 2>&1 || { echo "R module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v $cdo >/dev/null 2>&1 || { echo "CDO module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
$Rscript $PROGDIR/config/installpack.R


