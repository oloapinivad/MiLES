#!/bin/bash

# New pre-processor which can be used for any type of files
# It requires only that the files from the same variables are lying in the same folder
# Originally developed for Z500 only.

# Interpolation on regolar 2.5x2.5 grid, NH selection, daily averages, level extraction

# options
OPTIND=1
while getopts "h:p:d:e:r:v:l:o:c:f:" OPT; do
    case "$OPT" in
    h|\?) echo "Usage: universal_prepare.sh -p <project> -d <dataset> -e <expid> -r <ensemble>  \
                -v <varname> -l <levelout>  \
		-o <outputfilename> -c <configfile> -f <doforcedata>"
          exit 0 ;;
    p)    project=$OPTARG ;;
    d)    dataset=$OPTARG ;;
    e)    expid=$OPTARG ;;
    r)    ens=$OPTARG ;;
    v)    varname=$OPTARG ;;
    l)    levelout=$OPTARG ;;
    o)	  outputfilename=$OPTARG ;;
    c) 	  config=$OPTARG ;;
    f)    doforcedata=$OPTARG ;;
    esac
done
shift $((OPTIND-1))
echo $levelout

resolution=r144x73
#resolution=r288x145

if [[ ! -f $outputfilename ]] || [[ $doforcedata == "true" ]] ; then
	
	#clean old file
	rm -f $outputfilename
	identifier=${varname}${level}
	DATADIR=$(dirname $outputfilename)

	# clean ol tempdir and create a new one
	rm -rf $DATADIR/tempdir_${dataset}_${expid}_${ens}_*
    	TEMPDIR=$DATADIR/tempdir_${dataset}_${expid}_${ens}_$RANDOM
    	mkdir -p $TEMPDIR

	echo "${identifier} data are missing (or you forcily required a rebuilt!)"
	echo "... full data preparation is performed"
    	
	# machine dependent script (to call folder locations!)
	. config/config_${config}.sh

	# step 1: do it for NetCDF
	if [ -z ${expected_input_name} ] ; then 
		echo $INDIR/*.nc
		$cdonc cat $INDIR/*.nc $TEMPDIR/fullfile.nc
		# if NetCDF do not exists, check for grib files
		if [ $? -ne 0 ] ; then
			echo "perhaps you are using grib files..."
			$cdonc cat $INDIR/*.grb $TEMPDIR/fullfile.nc
			# if no file is found, exit the function with an error but not break the loop
			if [ $? -ne 0 ] ; then
				echo -e "${RED}ERROR: no file founds for this experiment! Exiting...${NC}"
				rm -rf ${DATADIR}
				return 1
			fi
		fi
	# otherwise use expected input name
	else
		#create a single huge file: not efficient but universal
		$cdonc cat $INDIR/${expected_input_name} $TEMPDIR/fullfile.nc
	fi

	# introducing a function to autofine geopotential height level in the file (looking for 500hPa or 50000Pa!)
	# updated with "level_select" command for special cases
	varunit=$( $cdonc zaxisdes $TEMPDIR/fullfile.nc | grep units | cut -f2 -d'"') 
	if [[ -z $varunit ]] ; then 
		echo -e "${RED}WARNING: Unknown unit for vertical axis!!!${NC}"
		echo "Selecting the first available level..."
	        echo "Assuming is $identifier ...it may be wrong!!!"
		level_select="-sellevidx,1"
	else
		if [[ "$varunit" == "millibar" ]] || [[ "$varunit" == "hPa" ]] ; then levelnew=$levelout ; fi
		if [[ "$varunit" == "Pa" ]] ; then levelnew=$((levelout*100)) ; fi
		echo "Level is $levelnew $varunit" 
		level_select="-sellevel,$levelnew"
	fi

	# main operations: sellevel and daymean + setlevel, setname, interpolation resolution and NH selection
	# CDO flag for extrapolation is on to avoid missing values for low res grid
	export REMAP_EXTRAPOLATE=on
	$cdonc sellonlatbox,0,360,0,90 -remapbil,$resolution -setlevel,${levelout}  \
		-setname,$varname -setunit,Pa -daymean ${level_select} \
		$TEMPDIR/fullfile.nc $TEMPDIR/smallfile.nc

	# new check: are there any missing values in the file?
	# use the CDO info output, slow
	missing=$(cdo info $TEMPDIR/smallfile.nc | awk -v x=1 '$7 >= x' | awk '{print $7}' | tail -n +3)
	echo $missing

	if [[ ! -z $missing ]] ; then
		echo -e "${RED}WARNING: $missing grid points have been found!!!${NC}"
		echo "Using CDO fillmiss to fill them, be careful!"
		$cdonc fillmiss $TEMPDIR/smallfile.nc $TEMPDIR/smallfile2.nc
		mv $TEMPDIR/smallfile2.nc $TEMPDIR/smallfile.nc
	fi


	# in order to avoid issues, all data are forced to be geopotential height 
	# in case geopotential is identified (i.e. values too large for a Z500
	sanityvalue=$($cdonc outputint -fldmean -seltimestep,1 $TEMPDIR/smallfile.nc)
	echo $sanityvalue

	#sanity check: convert to geopotential height if values are too high (to be improved for generalization)
	if [[ $sanityvalue -gt 10000 ]] ; then 
		echo "Values too high: considering as geopotential, dividing by g"
		$cdonc divc,9.80665 $TEMPDIR/smallfile.nc $TEMPDIR/smallfile2.nc
		mv $TEMPDIR/smallfile2.nc $TEMPDIR/smallfile.nc
	else
		echo "Geopotential height identified."
	fi	

	#copy to final file with absolute time axis which impoves R readability
	$cdo4 -a copy $TEMPDIR/smallfile.nc $outputfilename

    	#check cleaning
    	rm -f $TEMPDIR/*.nc
    	rmdir $TEMPDIR

else
	echo "${identifier} NetCDF data seems there, avoid universal_prepare.sh"
fi



