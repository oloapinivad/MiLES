################################################
#--MidLatitude Evaluation System for Ec-Earth--#
#------------------MiLES v0.3------------------#
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
INDIR=file.path("/home/paolo/work/DATA",exp,"day/Z500")

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
config="sansone"

################################################
system2(paste0(". config/config_",config,".sh"))
################################################

################################################
#-------------Z500 extraction------------------#
################################################

#call program for Z500 files: this program takes all the files
#into the $INDIR folder and prepare them in the single month files needed by MiLES
#since it is thought to be universal it is pretty much inefficient: it may be worth
#to personalize the script to obtain significant speedup
system2("diag_scripts/aux/miles/z500_prepare.sh",c(exp,toString(year1),toString(year2)))

for (season in seasons)
system2("Rscript "$PROGDIR/script/block_fast.R" $exp $year1 $year2 $season
        time $Rscript "$PROGDIR/script/block_figures.R" $exp $year1 $year2 $season
done
