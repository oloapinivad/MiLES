#!/bin/bash

# script for consecutive calls of wrapper_miles on different datasets, experiments or ensembles

project_exp=CMIP5
varname=zg

DIR=/work/users/paolo/data/${project_exp}
#DIR=$ARCHIVE/work/${project_exp}

#loop_dataset="ACCESS1-0  CMCC-CMS       EC-EARTH   GFDL-ESM2G  HadGEM2-CC_LT  IPSL-CM5A-MR  MIROC-ESM       MPI-ESM-P"
#ACCESS1-3     CanAM4   CMCC-CESM      CNRM-CM5       FGOALS-g2  GFDL-ESM2M  HadGEM2-ES     IPSL-CM5B-LR  MIROC-ESM-CHEM  MRI-CGCM3 
#bcc-csm1-1    CanESM2  CMCC-CESM_LOW  CSIRO-Mk3-6-0  FGOALS-s2  HadGEM2-A   inmcm4         MIROC4h       MPI-ESM-LR      MRI-ESM1 
#bcc-csm1-1-m  CCSM4    CMCC-CM        EC-Earth       GFDL-CM3   HadGEM2-CC  IPSL-CM5A-LR   MIROC5        MPI-ESM-MR      NorESM1-M"

loop_dataset=$(ls $DIR) 
#loop_dataset="ACCESS1-0 GFDL-ESM2M IPSL-CM5A-MR IPSL-CM5A-LR IPSL-CM5B-LR MRI-CGCM3"
expid="historical"


# mandatory
if [[ ${expid} == "rcp85" ]] ; then
	#year1_exp=2061
	#year2_exp=2100
	year1_exp=2006
	year2_exp=2100
elif [[ ${expid} == "historical" ]] ; then
	year1_exp=1950
	year2_exp=2005
	#year1_exp=1961
        #year2_exp=2000
elif [[ ${expid} == "sresa2" ]] ; then
        #year1_exp=2061
        #year2_exp=2100
	year1_exp=2046
        year2_exp=2100
fi



seasons="DJF"

# part of the code you want to run. Possible options are:
#"eofs": eofs parts
#"block": blocking
#"regimes": regimes
#"meandering": meandering index
#"figures": figures for all selected parts
#"forcedata": recompute the fullfile
#"forceanl" : recompute the analysis of each part
options="block figures"

for dataset_exp in $loop_dataset ; do
	loop_expid=$(ls $DIR/${dataset_exp}) 
	for expid_exp in $loop_expid ; do
		if [[ $expid_exp == $expid ]] ; then
			loop_ens=$(ls $DIR/${dataset_exp}/${expid_exp})
			for ens_exp in $loop_ens ; do
				echo $dataset_exp $ens_exp $expid_exp 
				. wrapper_miles.sh namelist/loop_namelist.sh
			done
		fi
	done
done
