#!/bin/bash



./rvage age -i ./truH.bin -o simres -m sim --sim simres.truth.txt -t 4 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 2.35e-08
./rvage age -i ./truH.bin -o simres -m sim --sim simres.truth.txt -t 4 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 2.35e-08
./rvage age -i ./truH.bin -o simres -m sim --sim simres.truth.txt -t 4 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 2.35e-08

