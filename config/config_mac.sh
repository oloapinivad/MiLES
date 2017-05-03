#!/bin/bash

#program definition
cdo=/opt/local/bin/cdo
cdonc="$cdo -f nc"
#if you CDO is not equipped of NetCDF4 compression change options here
cdo4="$cdo -f nc4 -z zip"
convert=/opt/local/bin/convert

#program folder
export PROGDIR=/Users/paolo/Desktop/newMiles/newMiles
#data folder
export DATADIR=/Users/paolo/Desktop/newMiles/data


####################################
#no need to change below this line #
####################################

#Z500 folder
export ZDIR=$DATADIR/Z500/$exp
#TEMPDIR folder
export TEMPDIR=$DATADIR/tempdir/${exp}_$RANDOM
#Blocking files folder
export BLOCKDIR=$DATADIR/files/Block/$exp
#NAO files folder
export EOFDIR=$DATADIR/files/EOFs/$exp
#figures folder
export FIGDIRBLOCK=$DATADIR/figures/Block/$exp
export FIGDIREOFS=$DATADIR/figures/EOFs/$exp

#creating folders
mkdir -p $ZDIR $FIGDIRBLOCK $FIGDIREOFS $BLOCKDIR $EOFDIR $TEMPDIR
echo $TEMPDIR
echo $BLOCKDIR
echo $FIGDIR

#safety check
echo "Check if R has been loaded"
command -v Rscript >/dev/null 2>&1 || { echo "R/3.1 module is not loaded. Aborting." >&2; exit 1; }
echo "R found: proceeding..."

echo "Check if NetCDF  has been loaded"
command -v ncdump >/dev/null 2>&1 || { echo "NetCDF module is not loaded. Aborting." >&2; exit 1; }
echo "NetCDF found: starting...."

#R check for key packages
Rscript $PROGDIR/config/installpack.R


