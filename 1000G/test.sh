#!/bin/bash

cd data

chr="20"

../rvage_dev --vcf ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --map genetic_map_GRCh37_chr${chr}.txt -o _test.1000G.chr${chr} &

cd ..

