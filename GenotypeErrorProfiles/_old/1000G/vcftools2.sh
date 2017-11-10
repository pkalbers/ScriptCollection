#!/bin/bash



for FILE in `ls ALL.*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
do

	vcftools --gzvcf $FILE \
	--indv NA12877 --remove-indels --recode \
	--out NA12877_${FILE} &

	sleep 5

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 4 ))
	do
        sleep 10
        JOBLIST=($(jobs -p))
    done

done
wait

