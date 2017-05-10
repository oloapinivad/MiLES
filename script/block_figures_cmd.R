#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]
FIGDIR=args[5]
EXPDIR=args[6]
REFDIR=args[7]
cfg=args[8]

miles.blockfigures<-function(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,EXPDIR,REFDIR,cfg)
