#!/bin/bash




./rvage ibd -o test.errG -i ../errH.bin -k 2 25 --randomSites 10 -m hmm -b 5000 --Ne 7300 -t 6 --hmm ../HMM.initial.prob.txt ./emiss.txt

