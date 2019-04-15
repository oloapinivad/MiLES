#!/bin/bash
#bash script to perform safety checks on ERAI 1979 - 2017
# sort of travis test zero order

TMPDIR=/scratch/users/paolo/safety_miles/
DATADIR=/work/users/paolo/miles/files
season=DJF
rm -rf $TMPDIR
mkdir -p $TMPDIR
cp -r  $DATADIR/Block/ERAI/1979_2017/DJF/*.nc $TMPDIR

cd ..
./wrapper_miles.sh namelist/erai_safety_check.R 
cd -

cdo diff $DATADIR/Block/ERAI/1979_2017/DJF/BlockFull_ERAI_1979_2017_DJF.nc $TMPDIR/BlockFull_ERAI_1979_2017_DJF.nc




