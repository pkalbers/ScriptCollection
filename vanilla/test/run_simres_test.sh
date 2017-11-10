#!/bin/bash



./rvage age -i ../vanilla.bin -o test -m sim --sim simres.test.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 0 -b 5000 --mut 1e-08
./rvage age -i ../vanilla.bin -o test -m sim --sim simres.test.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 0 -b 5000 --mut 1e-08
./rvage age -i ../vanilla.bin -o test -m sim --sim simres.test.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 0 -b 5000 --mut 1e-08

