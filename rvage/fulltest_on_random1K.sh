#!/bin/bash


./rvage --vcf OutOfAfricaHapMap20.H.vcf --map genetic_map_GRCh37_chr20.txt -o truH
./rvage --vcf OutOfAfricaHapMap20.P.vcf --map genetic_map_GRCh37_chr20.txt -o truP

./rvage --vcf OutOfAfricaHapMap20.GenErr_1000G.H.vcf --map genetic_map_GRCh37_chr20.txt -o errH
./rvage --vcf OutOfAfricaHapMap20.GenErr_1000G.P.vcf --map genetic_map_GRCh37_chr20.txt -o errP


# truth known

./rvage age -m dhg -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m dhg -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m dhg -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m dhg -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m fgt -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m fgt -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m fgt -i truH.bin -o truth.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m hmm -i truH.bin -o truth.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i truH.bin -o truth.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m hmm -i truH.bin -o truth.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m hmm -i truH.bin -o truth.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1


# truth phased

./rvage age -m dhg -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m dhg -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m dhg -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m dhg -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m fgt -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m fgt -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m fgt -i truP.bin -o truth.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m hmm -i truP.bin -o truth.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i truP.bin -o truth.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m hmm -i truP.bin -o truth.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m hmm -i truP.bin -o truth.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1



# error known

./rvage age -m dhg -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m dhg -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m dhg -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m dhg -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m fgt -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m fgt -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m fgt -i errH.bin -o error.H -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m hmm -i errH.bin -o error.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i errH.bin -o error.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m hmm -i errH.bin -o error.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m hmm -i errH.bin -o error.H -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1


# error phased

./rvage age -m dhg -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m dhg -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m dhg -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m dhg -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m fgt -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m fgt -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m fgt -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m fgt -i errP.bin -o error.P -b 5000 --positions pos1K.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1

./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 1
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 1 --useRecClock 0
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useMutClock 0 --useRecClock 1
./rvage age -m hmm -i errP.bin -o error.P -b 5000 --positions pos1K.txt --hmm HMM.initial.prob.txt HMM.emission.prob.txt -t 28 --Ne 7300 --Mr 2.35e-8 --useHardBreaks 1



