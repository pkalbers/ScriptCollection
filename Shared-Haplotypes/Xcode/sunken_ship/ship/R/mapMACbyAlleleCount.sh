#!/bin/bash

for FILE in $(ls ALL.*_integrated_*.vcf.gz)
do 
	screen -dm -S $FILE bash -c "Rscript mapMACbyAlleleCount.R $FILE 1000"
done
