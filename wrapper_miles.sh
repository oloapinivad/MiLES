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

# exp identificator: it is important for the folder structure.
# if you have more than on runs or experiments of the same model use
# this variable to distinguish them
# set also years and seasons to analyze
dataset_exp="CMCC-CMS"
year1_exp=1950
year2_exp=1970

# data folder: all the geopotential height data should be here
# NB: this is a folder structure used in my local machine
INDIR_EXP=/home/paolo/work/DATA/CMIP5/${dataset_exp}/HIST/r1/day/Z500
if [ "${dataset_exp}" == NCEP ] || [ "${dataset_exp}" == ERA40 ] || [ "${dataset_exp}" == ERAINTERIM  ] || [ "${dataset_exp}" == MERRA  ] ; then
	INDIR_EXP=/home/paolo/work/DATA/${dataset_exp}/day/Z500
fi


# std_clim flag: this is used to choose which climatology compare with results
# or with a user specified one: standard climatology is ERAINTERIM 1979-2014
# if std_clim=1 ERAINTERIM 1979-2014 is used
# if std_clim=0 a MiLES-generated different climatology can be specified
std_clim=1

# only valid if std_clim=0
dataset_ref="ERAINTERIM"
year1_ref=1980
year2_ref=2000
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
map_projection="no"
#map_projection="azequalarea"

#config name: create your own config file for your machine.
config=sansone


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
        year2_ref=2014
        REFDIR=$PROGDIR/clim
        exps=$dataset_exp
else
        REFDIR=$FILESDIR
        exps=$(echo ${dataset_exp} ${dataset_ref})
fi


################################################
#-------------Computation----------------------#
################################################

# loop to produce data: on experiment and - if needed - reference
for exp in $exps ; do

	# select for experiment
	if [[ $exp == $dataset_exp ]] ; then
		year1=${year1_exp}; year2=${year2_exp}; INDIR=${INDIR_EXP}
	fi

	# select for reference
	if [[ $exp == $dataset_ref ]] ; then
	        year1=${year1_ref}; year2=${year2_ref}; INDIR=${INDIR_REF}
	fi

	#definition of the fullfile name
	ZDIR=$OUTPUTDIR/Z500/$exp
	mkdir -p $ZDIR
	z500filename=$ZDIR/Z500_${exp}_fullfile.nc
	echo $z500filename

	#fullfile prepare
	time . $PROGDIR/script/z500_prepare.sh $exp $year1 $year2 $INDIR $z500filename
	for season in $seasons ; do
		echo $season
		# EOFs
		time . $PROGDIR/script/eof_fast.sh $exp $year1 $year2 "$seasons" $tele $z500filename $FILESDIR
		# blocking
		time $Rscript "$PROGDIR/script/block_fast.R" $exp $year1 $year2 $season $z500filename $FILESDIR $PROGDIR 
		# regimes
		time $Rscript "$PROGDIR/script/regimes_fast.R" $exp $year1 $year2 $season $z500filename $FILESDIR $PROGDIR $nclusters
	done

done
################################################
#-----------------Figures----------------------#
################################################

for season in $seasons ; do
	echo $season
	# EOFs figures
	time $Rscript "$PROGDIR/script/eof_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $tele
	# blocking figures
	time $Rscript "$PROGDIR/script/block_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR
	# regimes figures
	time $Rscript "$PROGDIR/script/regimes_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters
done






