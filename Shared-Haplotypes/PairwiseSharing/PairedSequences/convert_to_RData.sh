#!/bin/bash



for file in *.seqpair
do
	echo $file
	Rscript ./convert_to_RData.R $file &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 5 ))
	do
        sleep 10
        JOBLIST=($(jobs -p))
    done
done
wait

