#read command line
args <- commandArgs(TRUE)
exp=args[1]
year1=args[2]
year2=args[3]
season=args[4]
ZDIR=args[5]
BLOCKDIR=args[6]
PROGDIR=args[7]

source(paste0(PROGDIR,"/script/basis_functions.R"))
source(paste0(PROGDIR,"/script/block_fast.R"))
miles.blockfast(exp,year1,year2,season,ZDIR,BLOCKDIR)
