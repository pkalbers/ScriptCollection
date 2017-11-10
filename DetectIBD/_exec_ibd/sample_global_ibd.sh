#!/bin/bash

ROOT="../"
CORES=10

while true
do
	Rscript ${ROOT}_exec_msprime/sample_global_ibd.R &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $CORES )) ### max number of parallel processes
	do
        sleep $[ ( $RANDOM % 30 ) + 1 ]s
        JOBLIST=($(jobs -p))
    done
done

wait
