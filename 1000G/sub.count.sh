#!/bin/bash

for K in `seq 2 50`
do
	for C in `seq 1 22`
	do
		qsub -v chr=${C},frq=${K} ./job.count.sh
	done
done

