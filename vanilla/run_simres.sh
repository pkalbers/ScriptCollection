#!/bin/bash



#./rvage age -i ./vanilla.bin -o simres -m sim --sim simres.dis1000.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 1e-08
./rvage age -i ./vanilla.bin -o simres2 -m sim --sim simres.dis1000.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08
./rvage age -i ./vanilla.bin -o simres2 -m sim --sim simres.dis1000.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08

