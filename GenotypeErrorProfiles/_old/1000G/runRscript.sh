#!/bin/bash


for CHR in `seq 1 22`
do

	Rscript profileFromVcf.R \
	./platinum/NA12878.vcf.gz \
	./NA12878_ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.recode.vcf \
	./platinum/ConfidentRegions.bed \
	./ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.frq \
	profile.chr${CHR}.RData &

	sleep 5

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 4 ))
	do
        sleep 10
        JOBLIST=($(jobs -p))
    done

done
wait
