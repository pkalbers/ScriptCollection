#!/bin/bash

cd rvage

for PACK in `ls ../expos/`
do

echo
echo $PACK
echo

../rvage_dev age -i ../../vanilla.bin -o rvage_NN_${PACK} -m HMM --positions ../expos/$PACK -t 4 --Ne 10000 --mut 1e-08 --hmm ../_initials.hhmm.tru.txt ../_emission.hhmm.tru.100.txt --selectNN 1
../rvage_dev age -i ../../vanilla.bin -o rvage_RD_${PACK} -m HMM --positions ../expos/$PACK -t 4 --Ne 10000 --mut 1e-08 --hmm ../_initials.hhmm.tru.txt ../_emission.hhmm.tru.100.txt --selectNN 0

done
