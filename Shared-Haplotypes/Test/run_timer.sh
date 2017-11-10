#!/bin/bash


for i in `seq 2 25`;
do

	./ship_timer \
	-i ../ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
	--exact --threads 24 \
	-o exact.${i} -f ${i}

	./ship_timer \
	-i ../ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
	--threads 24 \
	-o cum.${i} -f ${i}

done



