#!/bin/bash

cd data

for chr in `seq 22 1`
do
	../rvage --vcf ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --map genetic_map_GRCh37_chr${chr}.txt -o 1000G.chr${chr} &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 4 ))
	do
        sleep 30
        JOBLIST=($(jobs -p))
    done
done

wait

