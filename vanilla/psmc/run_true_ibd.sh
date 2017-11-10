#!/bin/bash

cd true_ibd

for CCF in `ls ../rvage/rvage_*.ccf.txt | sort`
do

echo
echo $CCF
echo

Rscript ../prepare.true_ibd.R ${CCF}

FILE=`basename ${CCF}`

/usr/local/var/homebrew/linked/python3/bin/python3.6 ../msprime.true_ibd.py _${FILE} ../../vanilla.hdf5 > true_${FILE}

rm _${FILE}

done
