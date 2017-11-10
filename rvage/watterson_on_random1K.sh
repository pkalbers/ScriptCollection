#!/bin/bash


./rvage age -m dhg -i truH.bin -o WE.truth.H -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i truH.bin -o WE.truth.H -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i truH.bin -o WE.truth.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1

./rvage age -m dhg -i truP.bin -o WE.truth.P -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i truP.bin -o WE.truth.P -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i truP.bin -o WE.truth.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1

./rvage age -m dhg -i errH.bin -o WE.error.H -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i errH.bin -o WE.error.H -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i errH.bin -o WE.error.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1

./rvage age -m dhg -i errP.bin -o WE.error.P -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i errP.bin -o WE.error.P -b 5000 --positions pos1K.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i errP.bin -o WE.error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 4 --Ne 7300 --estTheta --useMutClock 1 --useRecClock 1

