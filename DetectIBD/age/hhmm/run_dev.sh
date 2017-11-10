#!/bin/bash



# ./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
#
# ./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
#
# ./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
#
# ./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
# ./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
#
# wait


./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 15 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 15 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 15 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 15 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 &
wait


./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./truH.bin -o result_truH -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &

./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./truP.bin -o result_truP -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &

./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./errH.bin -o result_errH -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &

./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 5 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 5 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &
./rvage_dev age -i ./errP.bin -o result_errP -m HMM --positions expos.equal.txt -t 6 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 500 -b 50000 --mut 2.35e-08 --hmmi &

wait
