#!/bin/bash
#$ -N UKB_hgmd
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _hgmd_err.txt
#$ -o _hgmd_log.txt
#$ -cwd


cd ./hgmd/


REF=${ref}
CHR=${chr}
POS=${pos}

IDX="hgmd.${REF}.chr${CHR}_${POS}.sharing.idx"


if [ ! -f "${IDX}" ]
then
    ../rvage/rvage -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin --position ${POS} --saveShareIndex
fi


../rvage/rvage age -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr
../rvage/rvage age -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

../rvage/rvage ibd -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m HMM --hmm ../rvage/HMM.initial.prob.txt ../rvage/HMM.emission.prob.txt --Ne 10000


../rvage/rvage age -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m DHG --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr
../rvage/rvage age -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m DHG --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

../rvage/rvage ibd -o hgmd.${REF}.chr${CHR}_${POS} -i ../data/ukb_chr${CHR}_c.bin -x ${IDX} -b 2000 -m DHG --writeIBD



