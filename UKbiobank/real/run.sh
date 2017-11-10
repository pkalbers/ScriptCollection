#!/bin/bash

REF="BRCA1"
CHR=17
POS=41247941



./rvage age -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin --position ${POS} -t 6 -b 2000 -m HMM --hmm HMM.initial.prob.txt HMM.emission.prob.txt --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr --saveShareIndex

IDX="${REF}.chr${CHR}_${POS}.sharing.idx"

#./rvage age -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m HMM --hmm HMM.initial.prob.txt HMM.emission.prob.txt --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

./rvage ibd -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin -x ${IDX} -t 6 -b 2000 -m HMM --hmm HMM.initial.prob.txt HMM.emission.prob.txt --Ne 10000


#./rvage age -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --useHardBreaks 0 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr
#./rvage age -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --useHardBreaks 1 --useMutClock 0 --useRecClock 1 --Ne 10000 --mut 1.2e-8 --writeDistr

#./rvage ibd -o ${REF}.chr${CHR}_${POS} -i ukb_chr${CHR}_c.bin -x ${IDX} -t 4 -b 2000 -m DHG --writeIBD







