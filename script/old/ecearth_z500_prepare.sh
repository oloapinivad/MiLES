#!/bin/bash

#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.
#parallelized on months

#define experiment and years
expname=$1
yy1=$2
yy2=$3

#replace experiment name in directories
INDIR2=$( echo ${INDIR} | sed -e "s|<EXP>|${expname}|g" )
ZDIR2=$( echo ${ZDIR} | sed -e "s|<EXP>|${expname}|g" )

#creating folder
mkdir -p $ZDIR2

for (( yy=$yy1; yy<=$yy2; yy++ )); do
	
	DDIR=$INDIR2/Output_${yy}/IFS
	echo $yy

	for mpar in $(seq 1 $nprocs 12) ; do
		for mm in $(seq $mpar $((mpar+nprocs-1)) ) ; do

		yyyymm=$( printf "%4d%02d" ${yy} ${mm} )
		if [ ! -f ${ZDIR2}/Z500_${expname}_${yyyymm}.nc ] || [ ${force_Z500} == 1 ] ; then

			 $cdo -f nc$nc_flag sellonlatbox,0,360,0,90 -remapcon2,r144x73 \
			-sp2gpl -daymean -shifttime,-3hours -sellevel,50000 \
                        -selcode,129 ${DDIR}/ICMSH$expname+${yyyymm} ${ZDIR2}/Z500_${expname}_${yyyymm}.nc &

		else
			#check if some is month is partially extracted
			nt=$( cdo -s ntime  ${ZDIR2}/Z500_${expname}_${yyyymm}.nc )
			if [ $nt -lt "28" ] ; then

				$cdo -f nc$nc_flag  sellonlatbox,0,360,0,90 -remapcon2,r144x73 \
				 -sp2gpl -daymean -shifttime,-3hours -sellevel,50000 \
				-selcode,129 ${DDIR}/ICMSH$expname+${yyyymm} ${ZDIR2}/Z500_${expname}_${yyyymm}.nc &
			fi
		fi
		done
		wait
	done
done


