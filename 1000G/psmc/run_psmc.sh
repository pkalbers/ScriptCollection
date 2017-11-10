#!/bin/bash

cd psmc

for CCF in `ls ../rvage/rvage_*.ccf.txt | sort`
do

echo
echo $CCF
echo

Rscript ../run_psmc.R ../rvage/${CCF}

done
