#!/bin/bash



#./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 1000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt
./rvage_dev age -i ../data/1000G.chr20.bin -o result_t1000 -m HMM --positions expos.chr20.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.1000.txt
./rvage_dev age -i ../data/1000G.chr20.bin -o result_t10000 -m HMM --positions expos.chr20.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.10000.txt


#
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08 --hmmi
#
#
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08 --hmmi
#
#
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08
# ./rvage_dev age -i ../data/1000G.chr20.bin -o result -m HMM --positions expos.chr20.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 1.2e-08 --hmmi
#

