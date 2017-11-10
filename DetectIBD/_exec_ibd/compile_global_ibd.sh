#!/bin/bash

ROOT="./"
PREFIX="OutOfAfricaHapMap20"
CORES=10

for IBDFILE in ${ROOT}global_ibd/global_ibd.*.RData
do
	echo ""
	echo $PREFIX
	echo $IBDFILE
	echo ""

	FLAG=`basename ${IBDFILE}`
	FLAG="flag.${FLAG}.txt"

	if [ -e $FLAG ]
	then
		echo "Already processed"
		continue
	fi

	> $FLAG

	Rscript ${ROOT}_exec_prepare/compile_global_ibd.R $PREFIX $IBDFILE &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $CORES )) ### max number of parallel processes
	do
        sleep $[ ( $RANDOM % 30 ) + 1 ]s
        JOBLIST=($(jobs -p))
    done
done

wait
