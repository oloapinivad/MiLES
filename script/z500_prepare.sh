#!/bin/bash

#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.

#define experiment and years
exp=$1
year1=$2
year2=$3
INDIR=$4
DATADIR=$5

TEMPDIR=$DATADIR/tempdir_${exp}_$RANDOM
ZDIR=$DATADIR/$exp
mkdir -p $TEMPDIR $ZDIR

#flag to force the analysis (0 let the code decide, 1 force analysis)
do_all=0

# check for existance of 12 months of data
for (( year=$year1; year<=$year2; year++ )); do
                nmonths=$(echo $(ls $ZDIR/Z500_${exp}_${year}*.nc | wc -l ))
                #echo $nmonths
                if [[ $nmonths != 12 ]] ; then do_all=1 ;  break ; fi
done

#check if needs to do all the analysis
if [[ ${do_all} -eq 1 ]] ; then

	echo "Z500 data are missing... full extraction is performed"

	#create a single huge file: not efficient but universal
	$cdonc cat $INDIR/*.nc $TEMPDIR/fullfile.nc
	$cdonc sellonlatbox,0,360,0,90 -remapcon2,r144x73 -setlevel,50000 -setname,zg -selyear,$year1/$year2 $TEMPDIR/fullfile.nc $TEMPDIR/smallfile.nc

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
	for (( year=$year1; year<=$year2; year++ )); do
			#yearyearmm=$( printf "%4d%02d" ${year} ${mm} )
			$cdo4 splitmon $TEMPDIR/Z500_year_${year}.nc $ZDIR/Z500_${exp}_${year}
	done

else
	#print that script is not run
	echo "All Z500 NetCDF data seems there, avoid z500_prepare.sh"
fi

#check cleaning
rm -f $TEMPDIR/*.nc
rmdir $TEMPDIR


