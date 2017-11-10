#!/bin/bash




./irvage age -i ./truH.bin -o tru_hap -m FGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m FGT --positions expos.10000.txt -t 60 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m FGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08


./irvage age -i ./truP.bin -o tru_phs -m FGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truP.bin -o tru_phs -m FGT --positions expos.10000.txt -t 60 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truP.bin -o tru_phs -m FGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08



./irvage age -i ./truH.bin -o tru_hap -m DGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m DGT --positions expos.10000.txt -t 60 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m DGT --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08


./irvage age -i ./truH.bin -o tru_hap -m HMM --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m HMM --positions expos.10000.txt -t 60 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
./irvage age -i ./truH.bin -o tru_hap -m HMM --positions expos.10000.txt -t 60 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08

