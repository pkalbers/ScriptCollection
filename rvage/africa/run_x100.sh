#!/bin/bash

TH=48
RS=100
RN=$RANDOM

METHODS="fgt|dhg|hmm --hmm HMM.initial.prob.txt HMM.emission.prob.txt"
DATASET="truH errH"

for DATA in `echo ${DATASET}`
do
	for MT in `seq 1 3`
	do
		METHOD=`echo ${METHODS} | cut -d'|' -f ${MT}`
		./rvage age -i ${DATA}.bin -k 2 26 --randomSites ${RS} -o x100.${DATA}.${RN} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --saveShareIndex
		./rvage age -i ${DATA}.bin -k 2 26 --randomSites ${RS} -o x100.${DATA}.${RN} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --saveShareIndex
		./rvage age -i ${DATA}.bin -k 2 26 --randomSites ${RS} -o x100.${DATA}.${RN} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --saveShareIndex
	done
done

