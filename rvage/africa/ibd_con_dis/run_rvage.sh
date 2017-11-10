#!/bin/bash


DATA="tru err"
PHSD="H P"

for DT in ${DATA}
do
for PH in ${PHSD}
do

FILE="ibd_con_dis.${DT}_${PH}.n5000_id"
./rvage_dev age -i ../${DT}${PH}.bin -o ${FILE}_NN -m HMM --positions n5000_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.${DT}.txt _emission.hhmm.${DT}.100.txt --selectNN 1
./rvage_dev age -i ../${DT}${PH}.bin -o ${FILE}_RD -m HMM --positions n5000_expos_id.txt -t 6 --Ne 7300 --mut 2.35e-08 --hmm _initials.hhmm.${DT}.txt _emission.hhmm.${DT}.100.txt --selectNN 0

done
done
