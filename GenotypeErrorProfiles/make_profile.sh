#!/bin/bash


cd ./1000g
Rscript ../prepare.error.R NA12878.chrALL.vcf ../_platinum/truth.NA12878.RData 1000g.weak.NA12878 > log.1000g.weak.NA12878.txt 2>&1 &
Rscript ../prepare.error.R NA12878.chrALL.vcf ../_platinum/truth.combined.NA12878.RData 1000g.strict.NA12878 > log.1000g.strict.NA12878.txt 2>&1 &
wait

cd ../affy6

Rscript ../prepare.error.R NA12878_ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz.recode.vcf ../_platinum/truth.NA12878.RData affy6.weak.NA12878 > log.affy6.weak.NA12878.txt 2>&1 &
Rscript ../prepare.error.R NA12878_ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz.recode.vcf ../_platinum/truth.combined.NA12878.RData affy6.strict.NA12878 > log.affy6.strict.NA12878.txt 2>&1 &

Rscript ../prepare.error.R NA12877_ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz.recode.vcf ../_platinum/truth.NA12877.RData affy6.weak.NA12877 > log.affy6.weak.NA12877.txt 2>&1 &
Rscript ../prepare.error.R NA12877_ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz.recode.vcf ../_platinum/truth.combined.NA12877.RData affy6.strict.NA12877 > log.affy6.strict.NA12877.txt 2>&1 &

wait


cd ../omni

Rscript ../prepare.error.R NA12878_ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz.recode.vcf ../_platinum/truth.NA12878.RData omni.weak.NA12878 > log.omni.weak.NA12878.txt 2>&1 &
Rscript ../prepare.error.R NA12878_ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz.recode.vcf ../_platinum/truth.combined.NA12878.RData omni.strict.NA12878 > log.omni.strict.NA12878.txt 2>&1 &

Rscript ../prepare.error.R NA12877_ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz.recode.vcf ../_platinum/truth.NA12877.RData omni.weak.NA12877 > log.omni.weak.NA12877.txt 2>&1 &
Rscript ../prepare.error.R NA12877_ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz.recode.vcf ../_platinum/truth.combined.NA12877.RData omni.strict.NA12877 > log.omni.strict.NA12877.txt 2>&1 &

wait

