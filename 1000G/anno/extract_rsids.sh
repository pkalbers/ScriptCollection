#!/bin/bash

for i in `seq 1 22`;
do
	gunzip -c ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz | grep -e 'rs[0-9]' > rsids.chr${i}.txt
done


