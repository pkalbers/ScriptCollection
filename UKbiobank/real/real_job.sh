#!/bin/bash
#$ -N UKB_real
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _real_err.txt
#$ -o _real_log.txt
#$ -cwd
#$ -pe shmem 4


cd ./real/




REF="---"
CHR=7
POS=107323983


../rvage/rvage age -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin --position ${POS} -t 4 -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr --saveShareIndex

IDX="real.chr{CHR}_${POS}.sharing.idx"

../rvage/rvage age -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

../rvage/rvage ibd -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --Ne 10000


../rvage/rvage age -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr
../rvage/rvage age -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

../rvage/rvage ibd -o real.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --writeIBD



