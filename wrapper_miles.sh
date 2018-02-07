#!/bin/bash
##
################################################
#--MidLatitude Evaluation System for Ec-Earth--#
#------------------MiLES v0.4------------------#
#----------Jun 2017, P. Davini, ISAC-CNR-------#
#
#
#
#
################################################
#- ------user configurations variables---------#
################################################

#config name: create your own config file for your machine.
config=sansone

# flag specific for ECMWF data structure
# it is used to call a specific preprocessing tool that follows 
# ensemble structure and folder at ECMWF
ECMWF=0

# exp identificator: it is important for the folder structure.
# if you have more than one runs (i.e. ensemble members) or experiments of the same model use
# this variable to distinguish them
# set also years 

#loop to create the ensembles
year1_exp=1979
year2_exp=2014
dataset_exp="CMCC-CM2"
#ens_list=$(seq -f "%02g" 24 24 )
ens_list="NO"

# INDIR_EXP ->data folder: all the geopotential height data should be here
# NB: this is a folder structure used in my local machine
# it is not used when you use ECMWF=1
#INDIR_EXP=/home/paolo/work/DATA/CMIP5/${dataset_exp}/HIST/r1/day/Z500
INDIR_EXP=/home/paolo/scratch/irene
if [ "${dataset_exp}" == NCEP ] || [ "${dataset_exp}" == ERA40 ] || [ "${dataset_exp}" == ERAINTERIM  ] || [ "${dataset_exp}" == MERRA  ] ; then
    INDIR_EXP=/home/paolo/work/DATA/${dataset_exp}/day/Z500
fi

# std_clim flag: this is used to choose which climatology compare with results
# or with a user specified one: standard climatology is ERAINTERIM 1979-2014
# if std_clim=1 ERAINTERIM 1979-2014 is used
# if std_clim=0 a MiLES-generated different climatology can be specified
std_clim=1

# only valid if std_clim=0
dataset_ref="S5"
ens_ref=mean
year1_ref=1982
year2_ref=2016
# NB: this is a folder structure used in my local machine
# it is not used when you use ECMWF=1
INDIR_REF=/home/paolo/work/DATA/CMIP5/${dataset_ref}/HIST/r1/day/Z500
if [ "${dataset_ref}" == NCEP ] || [ "${dataset_ref}" == ERA40 ] || [ "${dataset_ref}" == ERAINTERIM  ] || [ "${dataset_ref}" == MERRA  ] ; then
        INDIR_REF=/home/paolo/work/DATA/${dataset_ref}/day/Z500
fi

# please specify one or more of the 4 standard seasons using 3 characters
#seasons="DJF MAM SON JJA"
seasons="DJF"

# select which EOFs you want to compute
# "NAO": the 4 first  EOFs of North Atlantic, i.e. North Atlantic Oscillation as EOF1
# "AO" : the 4 first EOFs of Northern Hemispiere, i.e. Arctic Oscillation as EOF1 
# "lon1_lon2_lat1_lat2" : custom regions for EOFs: beware that std_clim will be set to 0!
tele="NAO"
#tele="-50_20_10_80"

# output file type for figures (pdf, png, eps)
# pdf are set by default
output_file_type="pdf"

# map projection that is used for plotting
# "no": standard lon-lat plotting (fastest)
# "azequalarea": polar plot with equal area
# these are suggested: any other polar plot by "mapproj" R package are supported
#map_projection="no"
map_projection="azequalarea"



###############################################
#-------------Configuration scripts------------#
################################################

# machine dependent script (set above)
. config/config_${config}.sh

# this script controls some of the graphical parameters
# as plot resolutions and palettes
CFGSCRIPT=$PROGDIR/config/config.R

# select how many clusters for k-means over the North Atlantic
# NB: only 4 clusters supported so far.  
nclusters=4


################################################
##NO NEED TO TOUCHE BELOW THIS LINE#############
################################################

#check if you can compare run with std_clim=1
if  [ $std_clim -eq 1 ] ; then
        if ! { [ "$tele" = NAO ] || [ "$tele" = AO ]; } ; then
                echo "Error: you cannot use non-standard EOFs region with std_clim=1"
                exit
        fi
fi

# if we are using standard climatology
if [[ ${std_clim} -eq 1 ]] ; then
    dataset_ref="ERAINTERIM"
    year1_ref=1979
    year2_ref=2016
    REFDIR=$PROGDIR/clim
    exps=$dataset_exp
	ens_ref="NO"
else
        REFDIR=$FILESDIR
        exps=$(echo ${dataset_exp} ${dataset_ref})
fi


################################################
#-------------Computation----------------------#
################################################

#ensemble loop
for ens_exp in ${ens_list} ; do

echo ${ens_exp}
echo $exps

# loop to produce data: on experiment and - if needed - reference
for exp in $exps ; do


	# select for experiment
	if [[ $exp == $dataset_exp ]] ; then
		year1=${year1_exp}; year2=${year2_exp}; INDIR=${INDIR_EXP}; ens=${ens_exp}
	fi

	# select for reference
	if [[ $exp == $dataset_ref ]] ; then
		if [ "${ens_ref}" == "mean" ] ; then continue ; fi
		if [ ${std_clim} -eq 1 ] ; then
 			echo "skip!"
		else
	        	year1=${year1_ref}; year2=${year2_ref}; INDIR=${INDIR_REF}; ens=${ens_ref}
		fi
	fi
	
	echo $exp $year1 $year2

	#definition of the fullfile name
	ZDIR=$OUTPUTDIR/Z500/$exp
	mkdir -p $ZDIR
	if [ "$ens" == "NO" ] ; then
		z500filename=$ZDIR/Z500_${exp}_fullfile.nc
	else
		z500filename=$ZDIR/${ens}/Z500_${exp}_${ens}_fullfile.nc
	fi
	echo $z500filename

	#fullfile prepare
	if [ $ECMWF -eq 0 ] ; then
		time . $PROGDIR/script/z500_prepare.sh $exp $year1 $year2 $INDIR $z500filename
	else
		time . $PROGDIR/script/ecmwf_z500_prepare.sh $exp $ens $year1 $year2 $z500filename
	fi

	for season in $seasons ; do
		echo $season
		# EOFs
		time . $PROGDIR/script/eof_fast.sh $exp $ens $year1 $year2 "$seasons" $tele $z500filename $FILESDIR
		# blocking
		time $Rscript "$PROGDIR/script/block_fast.R" $exp $ens $year1 $year2 $season $z500filename $FILESDIR $PROGDIR 
		# regimes
		time $Rscript "$PROGDIR/script/regimes_fast.R" $exp $ens $year1 $year2 $season $z500filename $FILESDIR $PROGDIR $nclusters
	done

done
################################################
#-----------------Figures----------------------#
################################################

for season in $seasons ; do
	echo $season
	# EOFs figures
	time $Rscript "$PROGDIR/script/eof_figures.R" $dataset_exp $ens_exp $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $tele
	# blocking figures
	time $Rscript "$PROGDIR/script/block_figures.R" $dataset_exp $ens_exp $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR
	# regimes figures
	time $Rscript "$PROGDIR/script/regimes_figures.R" $dataset_exp $ens_exp $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters
done

done

################################################
#--------------Ensemble mean-------------------#
################################################

if [ "${ens_list}" != "NO" ] ; then

	echo $FILESDIR
	for season in $seasons ; do
		MEANBLOCKDIR=$FILESDIR/${dataset_exp}/mean/Block/${year1}_${year2}/$season
		rm -rf $MEANBLOCKDIR; mkdir -p $MEANBLOCKDIR
		$cdo cat $FILESDIR/${dataset_exp}/*/Block/${year1}_${year2}/$season/BlockClim*.nc $MEANBLOCKDIR/full.nc
		$cdo timmean $MEANBLOCKDIR/full.nc $MEANBLOCKDIR/BlockClim_${dataset_exp}_mean_${year1}_${year2}_${season}.nc
		time $Rscript "$PROGDIR/script/block_figures.R" $dataset_exp mean $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR		

	done
fi

rm -f $PROGDIR/Rplots.pdf #remove sporious file creation by R	

