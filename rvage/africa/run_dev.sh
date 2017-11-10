#!/bin/bash



./rvage_dev age -i ./truH.bin -o new_tru_H_100_NN_raw -m HMM --positions n50_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.tru.txt _emission.hhmm.tru.100.txt --selectNN 1 --filterProp 0 --filterAuto 0
./rvage_dev age -i ./truH.bin -o new_tru_H_100_NN_adj -m HMM --positions n50_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.tru.txt _emission.hhmm.tru.100.txt --selectNN 1 --filterProp 0 --filterAuto 1

./rvage_dev age -i ./errH.bin -o new_err_H_100_NN_raw -m HMM --positions n50_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 0
./rvage_dev age -i ./errH.bin -o new_err_H_100_NN_adj -m HMM --positions n50_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.100.txt --selectNN 1 --filterProp 0 --filterAuto 1



# TIME="100 1000 10000"
#
# for TM in ${TIME}
# do
#
# echo
# echo $TM
# echo
#
#
# ./rvage_dev age -i ./errH.bin -o nn_errH_id_t${TM} -m HMM --positions expos.f_id.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --mut 2.35e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.${TM}.txt
# ./rvage_dev age -i ./errH.bin -o nn_errH_fn_t${TM} -m HMM --positions expos.f_fn.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --mut 2.35e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.${TM}.txt
# ./rvage_dev age -i ./errH.bin -o nn_errH_fp_t${TM} -m HMM --positions expos.f_fp.txt -t 6 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --mut 2.35e-08 --hmm _initials.hhmm.err.txt _emission.hhmm.err.${TM}.txt
#
#
# done


#./rvage_dev age -i ./truH.bin -o dev.truH -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
#./rvage_dev age -i ./truH.bin -o dev.truH -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
#./rvage_dev age -i ./truH.bin -o dev.truH -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
#
#./rvage_dev age -i ./truP.bin -o dev.truP -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 1 --useRecClock 0 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
#./rvage_dev age -i ./truP.bin -o dev.truP -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 0 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08
#./rvage_dev age -i ./truP.bin -o dev.truP -m HMM --hmmi --writeIBD --positions expos.10000.txt -t 4 --useMutClock 1 --useRecClock 1 --useHardBreaks 1 --maxDiscordant 5000 -b 50000 --mut 2.35e-08


