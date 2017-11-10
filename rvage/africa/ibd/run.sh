#!/bin/bash


./rvage ibd -o all.truH -i ../truH.bin -k 2 25 -m fgt -b 1000 -t 4
./rvage ibd -o all.truP -i ../truP.bin -k 2 25 -m fgt -b 1000 -t 4
./rvage ibd -o all.truG -i ../truH.bin -k 2 25 -m dhg -b 1000 -t 4


./rvage ibd -o all.errH -i ../errH.bin -k 2 25 -m fgt -b 1000 -t 4
./rvage ibd -o all.errP -i ../errP.bin -k 2 25 -m fgt -b 1000 -t 4
./rvage ibd -o all.errG -i ../errH.bin -k 2 25 -m dhg -b 1000 -t 4


./rvage ibd -o all.truG -i ../truH.bin -k 2 25 -m hmm --hmm ../HMM.initial.prob.txt ../HMM.emission.prob.txt -b 1000 --Ne 7300 -t 6
./rvage ibd -o all.errG -i ../errH.bin -k 2 25 -m hmm --hmm ../HMM.initial.prob.txt ../HMM.emission.prob.txt -b 1000 --Ne 7300 -t 6

