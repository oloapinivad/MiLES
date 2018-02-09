#!/bin/bash

#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.

#define experiment and years
dataset=$1
ensemble=$2
year1=$3
year2=$4
z500filename=$5
config=$6

if [ ! -f $z500filename ] ; then
    
    DATADIR=$(dirname $z500filename)
    TEMPDIR=$DATADIR/tempdir_${dataset}_$RANDOM
    mkdir -p $TEMPDIR

	echo "Z500 data are missing... full data preparation is performed"
    
    # machine dependent script (to call folder locations!)
    . config/config_${config}.sh

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


