#!/bin/bash

for CHR in `seq 1 22`
do
	echo "Chromosome ${CHR}"
	FILE=`echo ALL.chr${CHR}.*.vcf.gz`
	echo $FILE
	./ship -i ./${FILE} -f 50 -s ./integrated_call_samples_v3.20130502.ALL.panel --threads 4 --remove_unknown_markers --remove_multiallelic -o chr${CHR}.f50
done

for CHR in `seq 1 22`
do
	echo "Chromosome ${CHR}"
	FILE=`echo ALL.chr${CHR}.*.vcf.gz`
	echo $FILE
	./ship -i ./${FILE} -f 100 -s ./integrated_call_samples_v3.20130502.ALL.panel --threads 4 --remove_unknown_markers --remove_multiallelic -o chr${CHR}.f100
done



