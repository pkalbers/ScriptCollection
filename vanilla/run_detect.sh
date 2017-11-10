#!/bin/bash



# ./rvage age -i ./vanilla.bin -o detect -m fgt --positions expos.10000.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
# ./rvage age -i ./vanilla.bin -o detect -m fgt --positions expos.10000.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
# ./rvage age -i ./vanilla.bin -o detect -m fgt --positions expos.10000.txt -t 5 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
#
#
# ./rvage age -i ./vanilla.bin -o detect -m dgt --positions expos.10000.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
# ./rvage age -i ./vanilla.bin -o detect -m dgt --positions expos.10000.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
# ./rvage age -i ./vanilla.bin -o detect -m dgt --positions expos.10000.txt -t 5 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000


./rvage age -i ./vanilla.bin -o detect -m hmm --positions expos.10000.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
./rvage age -i ./vanilla.bin -o detect -m hmm --positions expos.10000.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000
./rvage age -i ./vanilla.bin -o detect -m hmm --positions expos.10000.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 -b 5000 --mut 1e-08 --maxDiscordant 1000

