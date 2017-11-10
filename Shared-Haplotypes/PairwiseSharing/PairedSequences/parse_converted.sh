#!/bin/bash

for file in subsample.*.seqpair.RData
do
	echo $file
	Rscript ./parse_converted.R $file &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 5 ))
	do
        sleep 10
        JOBLIST=($(jobs -p))
    done
done
wait

