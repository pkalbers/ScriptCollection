#!/bin/bash

for file in *.seqpair
do
	echo $file
	Rscript ./read_and_parse.R $file &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 3 ))
	do
        sleep 10
        JOBLIST=($(jobs -p))
    done
done
wait

