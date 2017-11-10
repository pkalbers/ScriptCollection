#!/bin/bash

VCF_H=`ls *.H.vcf`
OUT_H=${VCF_H%.*}
OUT_H="beagle_ibd.${OUT_H}"

VCF_P=`ls *.P.vcf`
OUT_P=${VCF_P%.*}
OUT_P="beagle_ibd.${OUT_P}"

GMAP="plink.chr20.GRCh37.map"

echo "Running Beagle IBD"
echo "   ${VCF_H}"
echo "   ${VCF_P}"
echo
echo "Output file prefix"
echo "   ${OUT_H}"
echo "   ${OUT_P}"
echo
echo "Genetic Map: ${GMAP}"
echo

java -jar ./beagle.16Jun16.7e4.jar gt=${VCF_H} ibd=true out=${OUT_H} map=${GMAP} nthreads=6 > log.${OUT_H}.txt &
java -jar ./beagle.16Jun16.7e4.jar gt=${VCF_P} ibd=true out=${OUT_P} map=${GMAP} nthreads=6 > log.${OUT_P}.txt
wait

echo
echo "DONE"
echo
