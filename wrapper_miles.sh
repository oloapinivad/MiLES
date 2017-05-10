#!/bin/bash

################################################
#--MidLatitude Evaluation System for Ec-Earth--#
#------------------MiLES v0.2------------------#
#----------May 2017, P. Davini, ISAC-CNR-------#
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
exp="ERAINTERIM"

# data folder: all the geopotential height data should be here
#INDIR=/home/paolo/work/DATA/CMIP5/$exp/AMIP/r1/day/Z500
INDIR=/home/paolo/work/DATA/$exp/day/Z500

#years and seasons to analyze
year1=1979
year2=2014

#please specify one or more of the 4 standard seasons using 3 characters
#seasons="DJF MAM JJA SON"
seasons="DJF"

#select NAO/AO
teles="NAO"

#output file type for figures (pdf, png, eps)
output_file_type="pdf"

#config name: create your own config file for your machine.
config=sansone

################################################
. config/config_${config}.sh
################################################


################################################
#-------------Z500 extraction------------------#
################################################

#call program for Z500 files: this program takes all the files
#into the $INDIR folder and prepare them in the single month files needed by MiLES
#since it is thought to be universal it is pretty much inefficient: it may be worth
#to personalize the script to obtain significant speedup
#. $PROGDIR/script/z500_prepare.sh $exp $year1 $year2

################################################
#-------EOFs computation and figures-----------#
################################################

#call to program for EOFs index/pattern. CDO-based, fast and efficient
#NAO is the only one implemented for now
#figures are done using linear regressions of PCs on monthly anomalies
#against standard period for Reanalysis.

#time . $PROGDIR/script/eof_fast.sh $exp $year1 $year2 "$seasons" "$teles"
for tele in $teles ; do
	for season in $seasons ; do
		echo $season $tele
#		time $Rscript "$PROGDIR/script/eof_figures.R" $exp $year1 $year2 $season $tele
	done
done

################################################
#------Blocking Computation and Figures--------#
################################################

for season in $seasons ; do
	time $Rscript "$PROGDIR/script/block_fast.R" $exp $year1 $year2 $season
	time $Rscript "$PROGDIR/script/block_figures.R" $exp $year1 $year2 $season
done

################################################
# cleaning 
rm -f $TEMPDIR/*.nc
rm -r $TEMPDIR


