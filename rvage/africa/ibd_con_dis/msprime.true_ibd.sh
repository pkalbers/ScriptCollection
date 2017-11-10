#!/bin/bash


for INPUT in `ls __ibd_con_dis.tru*.ccf.txt`
do
echo $INPUT
OUTPUT="true${INPUT}"
/usr/local/Cellar/python3/3.6.1/bin/python3 msprime.true_ibd.py $INPUT /Users/pkalbers/Research/DetectIBD/OutOfAfricaHapMap20-v3.hdf5 > ${OUTPUT} &
done
wait

for INPUT in `ls __ibd_con_dis.err*.ccf.txt`
do
echo $INPUT
OUTPUT="true${INPUT}"
/usr/local/Cellar/python3/3.6.1/bin/python3 msprime.true_ibd.py $INPUT /Users/pkalbers/Research/DetectIBD/OutOfAfricaHapMap20-v3.hdf5 > ${OUTPUT} &
done
wait
