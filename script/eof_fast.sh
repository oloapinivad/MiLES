#!/bin/bash

expname=$1
yy1=$2
yy2=$3
seasons=$4
teles=$5

#preparing unique netcdf file
rm -f $EOFDIR/*.nc
rm -f $TEMPDIR/*.nc
$cdonc cat $ZDIR/Z500*nc $TEMPDIR/daily_file.nc
$cdonc monmean -selyear,$yy1/$yy2 $TEMPDIR/daily_file.nc $TEMPDIR/monthly_file.nc

for tele in $teles ; do
for season in $seasons ; do
	echo $season

	#fix folders and file names
	EOFDIR2=$EOFDIR/${tele}/${yy1}_${yy2}/${season}
	mkdir -p $EOFDIR2
        suffix=${expname}_${yy1}_${yy2}_${season}

	#select seasons, compute monthly anomalies
        $cdonc selseas,$season $TEMPDIR/monthly_file.nc $TEMPDIR/season_monthly.nc
        $cdonc -r ymonmean $TEMPDIR/season_monthly.nc $TEMPDIR/season_mean.nc
        $cdo4 -r -b 64 ymonsub $TEMPDIR/season_monthly.nc $TEMPDIR/season_mean.nc $EOFDIR2/Z500_monthly_anomalies_${suffix}.nc

	#fix borders for EOFs
        if [ "${tele}" = NAO ]; then
                box=-90,40,20,85
        fi
        if [ "${tele}" = AO ]; then
                box=0,360,20,85
        fi

	#select box for anomalies
	$cdonc sellonlatbox,${box}  $EOFDIR2/Z500_monthly_anomalies_${suffix}.nc $TEMPDIR/box_anomalies_monthly.nc	

	# -----------------------------------------------#
	# NB: due to change in CDO code after version 1.6.4
	# eigenvectors need to be normalized and area weighting must be specified for eofcoeffs
	# if you use a previous version please be aware that inconsistencies may arise
	# also, PCs are not standardized (this is done in the plotting tool)
	# -----------------------------------------------#

	# compute EOFs
        $cdonc -r eof,4 $TEMPDIR/box_anomalies_monthly.nc $EOFDIR2/${tele}_Z500_eigenvalues_${suffix}.nc $TEMPDIR/pattern.nc

	# normalize eigenvectors
	$cdonc enlarge,$TEMPDIR/pattern.nc -sqrt -fldsum -sqr $TEMPDIR/pattern.nc $TEMPDIR/factor.nc
	$cdonc div $TEMPDIR/pattern.nc $TEMPDIR/factor.nc $EOFDIR2/${tele}_Z500_pattern_${suffix}.nc

	# evalute grid area weights (implicitely used in eofs) and compute principal components
	$cdo gridweights $TEMPDIR/box_anomalies_monthly.nc $TEMPDIR/ww.nc
        $cdonc -r eofcoeff $EOFDIR2/${tele}_Z500_pattern_${suffix}.nc -mul $TEMPDIR/ww.nc $TEMPDIR/box_anomalies_monthly.nc $EOFDIR2/${tele}_monthly_timeseries_${suffix}_

	#-deprecated-#
        #export CDO_WEIGHT_MODE=off
        #$cdonc -r eof,4 $TEMPDIR/box_anomalies_monthly.nc $EOFDIR2/${tele}_Z500_eigenvalues_${suffix}.nc $TEMPDIR/pattern.nc
        #$cdonc enlarge,$TEMPDIR/pattern.nc -sqrt -fldsum -sqr $TEMPDIR/pattern.nc $TEMPDIR/factor.nc
        #$cdonc div $TEMPDIR/pattern.nc $TEMPDIR/factor.nc $EOFDIR2/${tele}_Z500_pattern_${suffix}.nc
        #$cdonc -r eofcoeff $EOFDIR2/${tele}_Z500_pattern_${suffix}.nc $TEMPDIR/box_anomalies_monthly.nc $EOFDIR2/${tele}_monthly_timeseries_${suffix}_ 
        #export CDO_WEIGHT_MODE=on


done
done
