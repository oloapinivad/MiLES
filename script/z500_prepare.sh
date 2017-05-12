#!/bin/bash

#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.

#define experiment and years
expname=$1
yy1=$2
yy2=$3
INDIR=$4
DATADIR=$5

TEMPDIR=$DATADIR/tempdir_${exp}_$RANDOM
DATADIR=$ZDIR/$exp
mkdir -p $TEMPDIR $ZDIR

#create a single huge file: not efficient but universal
$cdonc cat $INDIR/*.nc $TEMPDIR/fullfile.nc
$cdonc sellonlatbox,0,360,0,90 -remapcon2,r144x73 -setlevel,50000 -setname,zg -selyear,$yy1/$yy2 $TEMPDIR/fullfile.nc $TEMPDIR/smallfile.nc

#in order to avoid issues, all data are forced to be geopotential height in case geopotential is identified (i.e. values too large for a Z500
sanityvalue=$($cdonc outputint -fldmean -seltimestep,1 $TEMPDIR/smallfile.nc)
echo $sanityvalue

#sanity check
if [[ $sanityvalue -gt 10000 ]] ; then 
	echo "Values too high: considering as geopotential, dividing by g"
	$cdonc divc,9.80665 $TEMPDIR/smallfile.nc $TEMPDIR/smallfile2.nc
	mv $TEMPDIR/smallfile2.nc $TEMPDIR/smallfile.nc
else
	echo "Geopotential height identified."
fi

#splitting year and months
$cdonc splityear $TEMPDIR/smallfile.nc $TEMPDIR/Z500_year_
for (( yy=$yy1; yy<=$yy2; yy++ )); do
		#yyyymm=$( printf "%4d%02d" ${yy} ${mm} )
		$cdo4 splitmon $TEMPDIR/Z500_year_${yy}.nc ${DATADIR}/Z500_${expname}_${yy}
done

rm -f $TEMPDIR/*.nc
rmdir $TEMPDIR


