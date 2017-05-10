#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
dataset_ref=args[4]
year1_ref=args[5]
year2_ref=args[6]
season=args[7]
FIGDIR=args[8]
EXPDIR=args[9]
REFDIR=args[10]
cfg=args[11]
PROGDIR=args[12]

source(paste0(PROGDIR,"/script/basis_functions.R"))
source(paste0(PROGDIR,"/script/block_figures.R"))
miles.blockfigures(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,EXPDIR,REFDIR,cfg) 
