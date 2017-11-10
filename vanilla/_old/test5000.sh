#!/bin/bash


# ./rvage age -i vanilla.bin -o test5000 -t 6 -m dhg -b 1000 --positions ex.pos.5000.txt --useMutClock 1 --useRecClock 0 --useHardBreaks 0
./rvage age -i vanilla.bin -o test5000 -t 6 -m dhg -b 1000 --positions ex.pos.5000.txt --useMutClock 1 --useRecClock 1 --useHardBreaks 0
./rvage age -i vanilla.bin -o test5000 -t 6 -m dhg -b 1000 --positions ex.pos.5000.txt --useMutClock 1 --useRecClock 1 --useHardBreaks 1
./rvage age -i vanilla.bin -o test5000 -t 6 -m fgt -b 1000 --positions ex.pos.5000.txt --useMutClock 0 --useRecClock 1 --useHardBreaks 0

