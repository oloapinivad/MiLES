#!/bin/bash
set -eu
################################################
#------- MidLatitude Evaluation System --------#
#------------------MiLES v0.6------------------#
#----------May 2018, P. Davini, ISAC-CNR-------#
#
#
#
#
################################################
#- ------user configurations variables---------#
################################################

#config name: create your own config file for your machine.
machine=wilma

#control flags to check which sections should be run
doeof=false #EOFs section
doblock=false #Blocking section 
doregime=false #Regimes section
domeand=true #Meandering section
dofigs=true #Do you want figures?

#control flag for re-run of MiLES if files already exists
doforcedata=false
doforceanl=true

# dataset-experiment-ensemble identificators: important for the folder structure and the naming.
# if you have more than one run or experiments of the same model use this variable to distinguish them
# Beware that $ens_list is a list of ensemble on which you can run a loop
# set also years 
year1_exp=2039
year2_exp=2068
dataset_exp=SPHINX
expid_exp=T255
ens_list="lfb2 lfb3 lfb4"
#ens_list=NO
#ens_list=$(seq -f e"%02g" 1 4 )

# std_clim flag: this is used to choose which climatology compare with results
# or with a user specified one: standard climatology is ERAINTERIM 1979-2017
# if std_clim=true ERAINTERIM 1979-2017 is used
# if std_clim=false a MiLES-generated different climatology can be specified
std_clim=true

# only valid if std_clim=false
dataset_ref=MPI-ESM-P
expid_ref=HIST
ens_ref=r1
year1_ref=1951
year2_ref=2005

# please specify one or more of the 4 standard seasons using 3 characters. 
# std_clim is supported for these 4 seasons only. 
# To analyse the whole year use "ALL"
# Beta: now you can define your own season putting together 3-character string for each consecutive month you 
# want to include, for example "Jan_Feb_Mar".
seasons="DJF JJA"
#seasons="DJF"

# select which teleconnection pattern EOFs you want to compute
# "NAO": the 4 first  EOFs of North Atlantic, i.e. North Atlantic Oscillation as EOF1
# "AO" : the 4 first EOFs of Northern Hemispiere, i.e. Arctic Oscillation as EOF1 
# "PNA": the 4 first EOFs of North Pacific, i.e. Pacific North American Pattern as EOF1 (beta)
# "lon1_lon2_lat1_lat2" : custom regions for EOFs: beware that std_clim will be set to false!
#teles="NAO AO PNA"
teles="NAO"
#tele="-50_20_10_80"

# output file type for figures (pdf, png, eps)
# pdf are set by default
#output_file_type="pdf"

# map projection that is used for plotting
# "no": standard lon-lat plotting (fastest)
# "azequalarea": polar plot with equal area
# these are suggested: any other polar plot by "mapproj" R package are supported
# "azequalarea" set by default
#map_projection="azequalarea"

###############################################
#-------------Configuration scripts------------#
################################################

# machine dependent script (set above)
dataset=""; expid="" ens=""; expected_input_name="" #declared to avoid problems with set -u
. config/config_${machine}.sh

# this script controls some of the graphical parameters
# as plot resolutions and palettes
CFGSCRIPT=$PROGDIR/config/R_config.R

# select how many clusters for k-means over the North Atlantic
# NB: only 4 clusters supported so far.  
nclusters=4


################################################
##NO NEED TO TOUCHE BELOW THIS LINE#############
################################################

#check if you can compare run with std_clim=true
if ${std_clim} ; then

	for tele in $teles ; do
        	if ! { [ "$tele" = NAO ] || [ "$tele" = AO ] || [ "$tele" = PNA ]; } ; then
        	        echo "Error: you cannot use non-standard EOFs region with std_clim=true"
                	exit
        	fi
	done
	
	for season in $seasons ; do
        	if ! { [ "$season" = DJF ] || [ "$season" = MAM ] || [ "$season" = SON ] || [ "$season" = JJA ]; } ; then 
                	echo "Error: you cannot use non-standard seasons with std_clim=true"
                	exit
        	fi        
	done

fi

# if we are using standard climatology
if ${std_clim} ; then

	dataset_ref="ERAI_clim"
	year1_ref=1979
	year2_ref=2017
	REFDIR=$PROGDIR/clim
	datasets=$dataset_exp
	expid_ref=NO
	ens_ref=NO
else

        REFDIR=$FILESDIR
        datasets=$(echo ${dataset_exp} ${dataset_ref})
fi


################################################
#-------------Computation----------------------#
################################################

#ensemble loop
for ens_exp in ${ens_list} ; do

# loop to produce data: on experiment and - if needed - reference
for dataset in $datasets ; do

	# select for experiment
	if [[ $dataset == $dataset_exp ]] ; then
		year1=${year1_exp}; year2=${year2_exp}; expid=${expid_exp}; ens=${ens_exp}
	fi

	# select for reference
	if [[ $dataset == $dataset_ref ]] ; then
		if [ "${ens_ref}" == "mean" ] ; then continue ; fi
		if ${std_clim} ; then
 			echo "skip!"
		else
	        	year1=${year1_ref}; year2=${year2_ref}; expid=${expid_ref}; ens=${ens_ref}
		fi
	fi
	
	echo $dataset $expid $ens $year1 $year2

	#definition of the fullfile name
	ZDIR=$OUTPUTDIR/Z500
	mkdir -p $ZDIR
	zfile=Z500 
	for dcode in $dataset $expid $ens ; do
		if [[ "$dcode" != "NO" ]] ; then
			zfile=${zfile}_${dcode}
			ZDIR=$ZDIR/$dcode
		fi
	done
	z500filename=$ZDIR/${zfile}_fullfile.nc
	
	echo $z500filename

	#fullfile prepare
	time . $PROGDIR/script/z500_prepare.sh $dataset $expid $ens $year1 $year2 $z500filename $machine $doforcedata

	for season in $seasons ; do
		echo $season
        	
		# Rbased EOFs
		if $doeof ; then
			for tele in $teles ; do
        			time $Rscript "$PROGDIR/script/Rbased_eof_fast.R" 	$dataset $expid $ens $year1 $year2 $season \
											$tele $z500filename $FILESDIR $PROGDIR $doforceanl
			done
		fi
		# blocking
		if $doblock ; then
			time $Rscript "$PROGDIR/script/block_fast.R" 	$dataset $expid $ens $year1 $year2 $season \
									$z500filename $FILESDIR $PROGDIR $doforceanl
		fi

		# regimes
		if [[ $doregime == "true" ]] && [[ $season == DJF ]] ; then
			time $Rscript "$PROGDIR/script/regimes_fast.R"	$dataset $expid $ens $year1 $year2 $season \
									$z500filename $FILESDIR $PROGDIR $nclusters $doforceanl
		fi
		# meandering index
        	if $domeand ; then
        	    time $Rscript "$PROGDIR/script/meandering_fast.R" 	$dataset $expid $ens $year1 $year2 $season \
									$z500filename $FILESDIR $PROGDIR $doforceanl
        	fi
	done

done
################################################
#-----------------Figures----------------------#
################################################

if $dofigs ; then 
for season in $seasons ; do
	echo $season

	# EOFs figures
	if $doeof ; then
		for tele in $teles ; do
    			time $Rscript "$PROGDIR/script/Rbased_eof_figures.R" 	$dataset_exp $expid_exp $ens_exp $year1_exp $year2_exp \
										$dataset_ref $expid_ref $ens_ref $year1_ref $year2_ref \
										$season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $tele
		done
	fi

	# blocking figures
	if $doblock ; then
		time $Rscript "$PROGDIR/script/block_figures.R" 	$dataset_exp $expid_exp $ens_exp $year1_exp $year2_exp \
									$dataset_ref $expid_ref $ens_ref $year1_ref $year2_ref \
									$season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR
	fi

	# regimes figures
	if [[ $doregime == "true" ]] && [[ $season == DJF ]] ; then
		time $Rscript "$PROGDIR/script/regimes_figures.R" 	$dataset_exp $expid_exp $ens_exp $year1_exp $year2_exp \
									$dataset_ref $expid_ref $ens_ref $year1_ref $year2_ref \
									$season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters
	fi
done
fi

done

################################################
#--------------Ensemble mean-------------------#
################################################

#check how many ensemble you have
aens=($ens_list); nens=${!ens[@]}
echo "Ensemble members are $ens_list"

#For blocking only
if [[ "${ens_list}" == "NO" ]] || [[ ${nens} -eq 0 ]] ; then
	echo "Only one ensemble member, exiting..."
else
	echo "Create ensemble mean... "
	#create loop using flags
	kinds=${kinds:-}
	[[ $doblock == true ]] && kinds="$kinds Block"   
	#[[ $doeof == true ]] && kinds="$kinds EOF" 
	#[[ $doregime == true ]] && kinds="$kinds Regimes" 
	echo $kinds

	for kind in $kinds ; do

		for season in ${seasons} ; do
		
			#case for file names
			case $kind in 
				Block) 	 
					MEANFILE=$FILESDIR/${dataset_exp}/mean/${kind}/${year1}_${year2}/$season/BlockClim_${dataset_exp}_mean_${year1}_${year2}_${season}.nc
					INFILES='$FILESDIR/${dataset_exp}/*/${kind}/${year1}_${year2}/$season/BlockClim_${dataset_exp}_*.nc'
					ARGS="$dataset_exp mean $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR"
					SCRIPT=$PROGDIR/script/block_figures.R
				;;	
				#Regimes)
				#        MEANFILE=$FILESDIR/${dataset_exp}/mean/${kind}/${year1}_${year2}/$season/RegimesPattern_${dataset_exp}_mean_${year1}_${year2}_${season}.nc
				#        INFILES='$FILESDIR/${dataset_exp}/*/${kind}/${year1}_${year2}/$season/RegimesPattern_${dataset_exp}_*.nc'
				#        ARGS="$dataset_exp mean $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters"
				#        SCRIPT=$PROGDIR/script/regimes_figures.R
				#	;;
			esac
	
			#mean and plot
			MEANDIR=$(dirname $MEANFILE)
			rm -rf $MEANDIR; mkdir -p $MEANDIR
			$cdo4 timmean -cat $(eval echo $INFILES) $MEANFILE
			time $Rscript "$SCRIPT" $ARGS		
			
		done
	done
fi

rm -f $PROGDIR/Rplots.pdf #remove sporious file creation by R	

