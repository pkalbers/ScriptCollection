#!/bin/bash



./rvage age -i ./vanilla.bin -o simres -m sim --sim truth.disALL.age.pairs.FGT.Mut.HardBreaks.txt -t 4 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 1e-06
./rvage age -i ./vanilla.bin -o simres -m sim --sim truth.disALL.age.pairs.FGT.Mut.HardBreaks.txt -t 4 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-06
./rvage age -i ./vanilla.bin -o simres -m sim --sim truth.disALL.age.pairs.FGT.Mut.HardBreaks.txt -t 4 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-06

