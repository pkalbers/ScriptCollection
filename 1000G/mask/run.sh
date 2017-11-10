#!/bin/bash


for CHR in `seq 1 22`;
do

FILE="test_1kg_chr${CHR}"

./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_NN_none -m HMM --positions _n10.expos.chr${CHR}.txt -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 0
./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_NN_auto -m HMM --positions _n10.expos.chr${CHR}.txt -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 1

./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_RD_none -m HMM --positions _n10.expos.chr${CHR}.txt -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 0 --filterProp 0 --filterAuto 0
./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_RD_auto -m HMM --positions _n10.expos.chr${CHR}.txt -t 6 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 0 --filterProp 0 --filterAuto 1

done

CHR="2"

FILE="__run_1kg_chr${CHR}"

./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_NN_none -m HMM --positions expos.chr2.lac.txt -t 4 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 0 --maxConcordant 1000 --maxDiscordant 1000
./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_NN_auto -m HMM --positions expos.chr2.lac.txt -t 4 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 1 --maxConcordant 1000 --maxDiscordant 1000

./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_RD_none -m HMM --positions expos.chr2.lac.txt -t 4 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 0 --filterProp 0 --filterAuto 0 --maxConcordant 1000 --maxDiscordant 1000
./rvage_dev age -i ../data/1000G.chr${CHR}.bin -o ${FILE}_RD_auto -m HMM --positions expos.chr2.lac.txt -t 4 --Ne 10000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 0 --filterProp 0 --filterAuto 1 --maxConcordant 1000 --maxDiscordant 1000

