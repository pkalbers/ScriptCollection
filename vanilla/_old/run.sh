#!/bin/bash

./rvage age -i vanilla.bin -o test -t 6 -m fgt -b 1000 -k 2 25 --randomSites 1000 --useMutClock 1 --useRecClock 1
./rvage age -i vanilla.bin -o test -t 6 -m fgt -b 1000 -k 2 25 --randomSites 1000 --useMutClock 1 --useRecClock 0
./rvage age -i vanilla.bin -o test -t 6 -m fgt -b 1000 -k 2 25 --randomSites 1000 --useMutClock 0 --useRecClock 1
./rvage age -i vanilla.bin -o test -t 6 -m fgt -b 1000 -k 2 25 --randomSites 1000 --useHardBreaks 1

