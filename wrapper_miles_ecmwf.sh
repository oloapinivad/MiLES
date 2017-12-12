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

#flag specific for ECMWF data structure
ECMWF=1


#loop to create the ensembles
#for ens in $(seq -f "%02g" 0 0 ) ; do
#	dataset_list=$(echo $dataset_list S5_ens$ens)
#done 

dataset_list="S5_ens01"
for dataset_exp in ${dataset_list} ; do

year1_exp=1982
year2_exp=2016
#INDIR_EXP=/scratch/rd/nedd/regimes/ERAI_daily

# INDIR_EXP ->data folder: all the geopotential height data should be here
# NB: this is a folder structure used in my local machine
# std_clim flag: this is used to choose which climatology compare with results
# or with a user specified one: standard climatology is ERAINTERIM 1979-2014
# if std_clim=1 ERAINTERIM 1979-2014 is used
# if std_clim=0 a MiLES-generated different climatology can be specified
std_clim=0

# only valid if std_clim=0
dataset_ref="S4_ens01"
year1_ref=1982
year2_ref=2016
INDIR_REF=/scratch/rd/nedd/regimes/ERAI_daily

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

#config name: create your own config file for your machine.
config=ecmwf


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
	if [ $ECMWF -eq 0 ] ; then
		time . $PROGDIR/script/z500_prepare.sh $exp $year1 $year2 $INDIR $z500filename
	else
		time . $PROGDIR/script/ecmwf_z500_prepare.sh $exp $year1 $year2 $z500filename
	fi

	for season in $seasons ; do
		echo $season
		# EOFs
		#time . $PROGDIR/script/eof_fast.sh $exp $year1 $year2 "$seasons" $tele $z500filename $FILESDIR
		# blocking
		time $Rscript "$PROGDIR/script/block_fast.R" $exp $year1 $year2 $season $z500filename $FILESDIR $PROGDIR 
		# regimes
		#time $Rscript "$PROGDIR/script/regimes_fast.R" $exp $year1 $year2 $season $z500filename $FILESDIR $PROGDIR $nclusters
	done

done
################################################
#-----------------Figures----------------------#
################################################

for season in $seasons ; do
	echo $season
	# EOFs figures
	#time $Rscript "$PROGDIR/script/eof_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $tele
	# blocking figures
	time $Rscript "$PROGDIR/script/block_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR
	# regimes figures
	#time $Rscript "$PROGDIR/script/regimes_figures.R" $dataset_exp $year1_exp $year2_exp $dataset_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters
done


done



