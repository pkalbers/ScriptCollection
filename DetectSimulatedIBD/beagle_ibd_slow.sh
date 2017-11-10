#!/bin/bash

VCF_H=`ls *.H.vcf`
OUT_H=${VCF_H%.*}
OUT_H="beagle_ibd.${OUT_H}"

VCF_P=`ls *.P.vcf`
OUT_P=${VCF_P%.*}
OUT_P="beagle_ibd.${OUT_P}"

echo "Running Beagle IBD"
echo "   ${VCF_H}"
echo "   ${VCF_P}"
echo
echo "Output file prefix"
echo "   ${OUT_H}"
echo "   ${OUT_P}"

java -jar /home/pkalbers/DetectSimulatedIBD/beagle.22Apr16.1cf.jar gt=${VCF_H} ibd=true out=${OUT_H} nthreads=2 > log.${OUT_H}.txt &
java -jar /home/pkalbers/DetectSimulatedIBD/beagle.22Apr16.1cf.jar gt=${VCF_P} ibd=true out=${OUT_P} nthreads=2 > log.${OUT_P}.txt
wait

echo
echo "DONE"
echo
