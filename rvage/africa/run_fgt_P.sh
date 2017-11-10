#!/bin/bash

METHOD="fgt"
TH=6


DATA="truP"

./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1
./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0
./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1



DATA="errP"

./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 1
./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 1 --useRecClock 0
./rvage age -i ${DATA}.bin -x ${DATA}.sharing.idx -o t1000.${DATA} -m $METHOD -t $TH -b 5000 --Ne 7300 --mut 1.2e-8 --useMutClock 0 --useRecClock 1

