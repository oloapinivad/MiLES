#!/bin/bash

#special prepare file for ECMWF
#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.

#define experiment and years
exp=$1
ens=$2
year1=$3
year2=$4
z500filename=$5

#extract dataset and ensemble for expname
#dataset=$(echo $exp | cut -d_ -f1)
#ensemble=$(echo $exp | cut -d"s" -f2)
dataset=$exp
ensemble=$ens

echo $dataset
echo $ensemble

if [ "${dataset}" == "ERA5" ] ; then
    INDIR=/scratch/rd/nedd/regimes/ERA5_daily
    expected_input_name=*ERA5*${ensemble}.grb
fi

if [ "${dataset}" == "ERAI" ] ; then
    INDIR=/scratch/rd/nedd/regimes/ERAI_daily
    expected_input_name=*ERAI*.grb
fi

if [ "${dataset}" == "S3" ] ; then
    INDIR=/home/ms/it/ccpd/regimes/S3_24h_nov
    expected_input_name=*/*S3*_${ensemble}.grb
fi

if [ "${dataset}" == "S4" ] ; then
    INDIR=/home/ms/it/ccpd/regimes/S4_daily_nov
    expected_input_name=*/*S4*_${ensemble}.grb
fi

if [ "${dataset}" == "S5" ] ; then
    INDIR=/scratch/rd/nedd/regimes/S5_24h_nov
    expected_input_name=*/*S5*_${ensemble}.grb
fi


if [ "${dataset}" == "S5LR" ] ; then
    INDIR=/scratch/rd/nedd/regimes/S5_24h_nov_lowres
    expected_input_name=*/*S5*_${ensemble}.grb
fi

echo $INDIR

DATADIR=$(dirname $z500filename)
TEMPDIR=$DATADIR/tempdir_${exp}_${ens}_$RANDOM
mkdir -p $TEMPDIR

if [ ! -f $z500filename ] ; then

	echo "Z500 data are missing... full data preparation is performed"

	#create a single huge file: not efficient but universal
	$cdonc cat $INDIR/${expected_input_name} $TEMPDIR/fullfile.nc

	# step 1: do it for NetCDF
	#$cdonc cat $INDIR/*.nc $TEMPDIR/fullfile.nc
	
	# if NetCDF do not exists, check for grib files
	#if [ $? -ne 0 ] ; then
	#	echo "perhaps you are using grib files..."
	#	$cdonc cat $INDIR/*.grb $TEMPDIR/fullfile.nc

	#fi
	
	#main operation: setlevel, name, resolution and NH
	$cdonc sellonlatbox,0,360,0,90 -remapcon2,r144x73 -setlevel,50000 -setname,zg $TEMPDIR/fullfile.nc $TEMPDIR/smallfile.nc

	#in order to avoid issues, all data are forced to be geopotential height in case geopotential is identified (i.e. values too large for a Z500
	sanityvalue=$($cdonc outputint -fldmean -seltimestep,1 $TEMPDIR/smallfile.nc)
	echo $sanityvalue

	#sanity check: convert to geopotential height
	if [[ $sanityvalue -gt 10000 ]] ; then 
		echo "Values too high: considering as geopotential, dividing by g"
		$cdonc divc,9.80665 $TEMPDIR/smallfile.nc $TEMPDIR/smallfile2.nc
		mv $TEMPDIR/smallfile2.nc $TEMPDIR/smallfile.nc
	else
		echo "Geopotential height identified."
	fi	

	#copy to final file with absolute time axis
	$cdo4 -a copy $TEMPDIR/smallfile.nc $z500filename

else
	echo "Z500 NetCDF data seems there, avoid z500_prepare.sh"
fi

#check cleaning
rm -f $TEMPDIR/*.nc
rmdir $TEMPDIR


