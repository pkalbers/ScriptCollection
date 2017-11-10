#!/bin/bash

for DATA in `echo OutOfAfricaHapMap20.H OutOfAfricaHapMap20.GenErr_1000G.H OutOfAfricaHapMap20.P OutOfAfricaHapMap20.GenErr_1000G.P`
do

./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0

./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m dhg -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1


./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0

./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m fgt -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1


./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt

./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o test.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt


./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt

./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i ${DATA}.vcf.bin -x ${DATA}.sharing.idx -o nopp.${DATA} -m hmm -t 7 -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt

done


