#!/bin/bash


./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt

./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --hmm HMM.initial.prob.txt HMM.emission.prob.txt


./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 

./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m fgt -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 


./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 

./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 
./rvage age -i errP.bin -x errP.sharing.idx -o t1000eP -m dhg -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 


./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt

./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt
./rvage age -i errP.bin -x errP.sharing.idx -o t1000ePnopp -m hmm -t 4 -b 1000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --usePosteriors 0 --hmm HMM.initial.prob.txt HMM.emission.prob.txt


