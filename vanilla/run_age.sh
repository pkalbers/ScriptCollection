#!/bin/bash


#./rvage -i ./vanilla.bin -o samplesites --saveShareIndex -k 2 11 --randomSites 100


./rvage age -i ./vanilla.bin -o check -m FGT --positions expos.10000.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 100 -b 5000 --mut 1e-06
./rvage age -i ./vanilla.bin -o check -m FGT --positions expos.10000.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 100 -b 5000 --mut 1e-06
./rvage age -i ./vanilla.bin -o check -m FGT --positions expos.10000.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 100 -b 5000 --mut 1e-06

