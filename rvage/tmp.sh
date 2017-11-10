#!/bin/bash



# error phased

./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 26 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 26 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 26 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 26 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1




