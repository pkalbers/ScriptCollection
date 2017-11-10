#!/bin/bash


PFX="PRDM9_more"
CHR="5"
POS="pos.PRDM9.more.txt"


for ITER in `seq 1 3`;
do

FILE="${PFX}_1kg_chr${CHR}.${ITER}"

./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}.NN -m HMM --positions $POS -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1
./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}.RD -m HMM --positions $POS -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 0

done

# --maxConcordant 1000 --maxDiscordant 1000
