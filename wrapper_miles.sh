#!/bin/bash
set -eu
################################################
#------- MidLatitude Evaluation System --------#
#------------------MiLES v0.7------------------#
#----------Nov 2018, P. Davini, CNR-ISAC-------#
#
#
#
#

usage()
{
   echo "Usage:"
   echo "       ./wrapper_miles.sh namelist"
   echo " 	Namelists are configuration files introduced since MiLES v0.7"
   echo "	You can specificy everything you need, a template is available namelist/namelist.tmpl"
   echo
}

if [ $# -ne 1 ]; then
   usage
   exit 1
fi

namelist=$1

if [ ! -f $namelist ] ; then
	"The namelist does not exists!"
	exit
else
	. $namelist
fi

# zero order configuration checker function
function has_config()
{
        for option in $options ; do
                if [[ $option == $1 ]] ; then
                        return 0
                        exit
                fi
        done
        return 1

}

################################################
#- ------user configurations variables---------#
################################################

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
project="" ; dataset=""; expid="" ens=""; expected_input_name="" #declared to avoid problems with set -u
. config/config_${machine}.sh

# this script controls some of the graphical parameters
# as plot resolutions and palettes
CFGSCRIPT=$PROGDIR/config/R_config.R

# select how many clusters for k-means over the North Atlantic
# NB: only 4 clusters supported so far.  
nclusters=4

# set variable identifier
identifier=${varname}${level}


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

# check if we need to force some part of the code
has_config forcedata && doforcedata=true || doforcedata=false
has_config forceanl && doforceanl=true || doforceanl=false

# set unbound variables
expid_exp=${expid_exp:-}
ens_exp=${ens_exp:-}
project_exp=${project_exp:-}
expid_ref=${expid_ref:-}
ens_ref=${ens_ref:-}
project_ref=${project_ref:-}


# if we are using standard climatology
if ${std_clim} ; then

	dataset_ref="ERAI_clim"
	year1_ref=1979
	year2_ref=2017
	REFDIR=$PROGDIR/clim
	datasets=$dataset_exp
	project_ref=""
	expid_ref=""
	ens_ref=""

else

        REFDIR=$FILESDIR
        datasets=$(echo ${dataset_exp} ${dataset_ref})
fi


################################################
#-------------Computation----------------------#
################################################

# loop to produce data: on experiment and - if needed - reference
for dataset in $datasets ; do

	# select for experiment
	if [[ $dataset == $dataset_exp ]] ; then
		year1=${year1_exp}; year2=${year2_exp}; expid=${expid_exp}; ens=${ens_exp}; project=${project_exp}
	fi

	# select for reference
	if [[ $dataset == $dataset_ref ]] ; then
		if ${std_clim} ; then
 			echo "skip!"
		else
	        	year1=${year1_ref}; year2=${year2_ref}; expid=${expid_ref}; ens=${ens_ref}; project=${project_ref}
		fi
	fi

	echo "year1=${year1}; year2=${year2}; expid=${expid}; ens=${ens}; project=${project}"
	
	# some informations on screen
	RED='\033[0;31m'; GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'
	echo	"###########################################################"
	echo 	"MiLES is running on:"
	echo -e	"${BLUE}$dataset $expid $ens $year1 $year2${NC}"
	echo 	"Active options are:"
        echo -e	"${BLUE}$options${NC}"
	echo    "###########################################################"

	#definition of the fullfile name
	ZDIR=$OUTPUTDIR/data/${identifier}/${project}
	mkdir -p $ZDIR
	zfile=${identifier}
	for dcode in $dataset $expid $ens ; do
		if [[ "$dcode" != "NO" ]] ; then
			zfile=${zfile}_${dcode}
			ZDIR=$ZDIR/$dcode
		fi
	done
	fullfilename=$ZDIR/${zfile}_fullfile.nc
	
	#echo $fullfilename

	#fullfile prepare
	time . $PROGDIR/script/universal_prepare_beta.sh -d "$dataset" -e "$expid" -r "$ens" -v $varname -l $level \
							-o $fullfilename -c $machine -f $doforcedata
	# exit if the prepare has gone wrong
	if [ "$?" == 1  ] ; then
		return 1
	fi
	#time . $PROGDIR/script/z500_prepare.sh $dataset $expid $ens $year1 $year2 $fullfilename $machine $doforcedata

	for season in $seasons ; do
		echo -e "${GREEN}Working on $season season${NC}"
        	
		# Rbased EOFs
		if has_config eofs ; then
			for tele in $teles ; do
        			time $Rscript "$PROGDIR/script/eof_fast.R" 	"$project" "$dataset" "$expid" "$ens" $year1 $year2 $season \
										$tele $fullfilename $FILESDIR $PROGDIR $doforceanl
			done
		fi
		# blocking
		if has_config block ; then
			[[ $varname == "zg" ]] && blockscript="block_fast.R"
			[[ $varname == "ua" ]] && blockscript="u500_block.R"
			time $Rscript "$PROGDIR/script/$blockscript" 	"$project" "$dataset" "$expid" "$ens" $year1 $year2 $season \
									$fullfilename $FILESDIR $PROGDIR $doforceanl
		fi

		# regimes
		if has_config regimes && [[ $season == DJF ]] ; then
			time $Rscript "$PROGDIR/script/regimes_fast.R"	"$project" "$dataset" "$expid" "$ens" $year1 $year2 $season \
									$fullfilename $FILESDIR $PROGDIR $nclusters $doforceanl
		fi
		# meandering index
        	if has_config meandering ; then
        	    	time $Rscript "$PROGDIR/script/meandering_fast.R" 	"$project" "$dataset" "$expid" "$ens" $year1 $year2 $season \
										$fullfilename $FILESDIR $PROGDIR $doforceanl
        	fi
	done

done
################################################
#-----------------Figures----------------------#
################################################

if has_config figures ; then 
for season in $seasons ; do
	echo     "###########################################################"
	echo -e  "${GREEN}Producing figures for $season season ${NC}"
	echo     "###########################################################"

	# EOFs figures
	if has_config eofs ; then
		for tele in $teles ; do
    			time $Rscript "$PROGDIR/script/eof_figures.R" 	"$project_exp" "$dataset_exp" "$expid_exp" "$ens_exp" "$year1_exp" "$year2_exp" \
									"$project_ref" "$dataset_ref" "$expid_ref" "$ens_ref" "$year1_ref" "$year2_ref" \
									$season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $tele
		done
	fi

	# blocking figures
	if has_config block ; then
		time $Rscript "$PROGDIR/script/block_figures.R" 	"$project_exp" "$dataset_exp" "$expid_exp" "$ens_exp" "$year1_exp" "$year2_exp" \
									"$project_ref" "$dataset_ref" "$expid_ref" "$ens_ref" "$year1_ref" "$year2_ref" \
									$varname $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR
	fi

	# regimes figures
	if has_config regimes && [[ $season == DJF ]] ; then
		time $Rscript "$PROGDIR/script/regimes_figures.R" 	"$project_exp" "$dataset_exp" "$expid_exp" "$ens_exp" "$year1_exp" "$year2_exp" \
								 	"$project_ref" "$dataset_ref" "$expid_ref" "$ens_ref" "$year1_ref" "$year2_ref" \
									$season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters
	fi
done
fi

################################################
#--------------Ensemble mean-------------------#
################################################

#check how many ensemble you have
#aens=($ens_list); nens=${!ens[@]}
#echo "Ensemble members are $ens_list"
#
##For blocking only
#if [[ "${ens_list}" == "NO" ]] || [[ ${nens} -eq 0 ]] ; then
	#echo "Only one ensemble member, exiting..."
#else
	#echo "Create ensemble mean... "
	##create loop using flags
	#kinds=${kinds:-}
	#[[ has_config block ]] && kinds="$kinds Block"   
	##[[ has_config eofs ]] && kinds="$kinds EOF" 
	##[[ has_config regimes ]] && kinds="$kinds Regimes" 
	#echo $kinds
#
	#for kind in $kinds ; do
#
		#for season in ${seasons} ; do
		#
			#case for file names
			#case $kind in 
				#Block) 	 
					#MEANFILE=$FILESDIR/${dataset_exp}/mean/${kind}/${year1}_${year2}/$season/BlockClim_${dataset_exp}_mean_${year1}_${year2}_${season}.nc
					#INFILES='$FILESDIR/${dataset_exp}/*/${kind}/${year1}_${year2}/$season/BlockClim_${dataset_exp}_*.nc'
					#ARGS="$dataset_exp mean $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR"
					#SCRIPT=$PROGDIR/script/block_figures.R
				#;;	
				#Regimes)
				##        MEANFILE=$FILESDIR/${dataset_exp}/mean/${kind}/${year1}_${year2}/$season/RegimesPattern_${dataset_exp}_mean_${year1}_${year2}_${season}.nc
				##        INFILES='$FILESDIR/${dataset_exp}/*/${kind}/${year1}_${year2}/$season/RegimesPattern_${dataset_exp}_*.nc'
				##        ARGS="$dataset_exp mean $year1_exp $year2_exp $dataset_ref $ens_ref $year1_ref $year2_ref $season $FIGDIR $FILESDIR $REFDIR $CFGSCRIPT $PROGDIR $nclusters"
				##        SCRIPT=$PROGDIR/script/regimes_figures.R
				##	;;
			#esac
	#
			##mean and plot
			#MEANDIR=$(dirname $MEANFILE)
			#rm -rf $MEANDIR; mkdir -p $MEANDIR
			#$cdo4 timmean -cat $(eval echo $INFILES) $MEANFILE
			#time $Rscript "$SCRIPT" $ARGS		
			#
		#done
	#done
#fi

rm -f $PROGDIR/Rplots.pdf #remove sporious file creation by R	

