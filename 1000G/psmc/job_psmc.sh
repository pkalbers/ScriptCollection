#!/bin/bash
#$ -N psmc_c20
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _psmc_c20_err.txt
#$ -o _psmc_c20_log.txt
#$ -cwd

cd psmc

for CCF in `ls ../rvage/rvage_*.ccf.txt | sort -R`
do

echo
echo $CCF
echo

FLAG=__flag__`basename ${CCF}`

sleep $[ ( $RANDOM % 10 )  + 1 ]s

if [ -f ${FLAG} ]
then
    continue
fi

sleep $[ ( $RANDOM % 10 )  + 1 ]s

if [ -f ${FLAG} ]
then
    continue
fi

> ${FLAG}

/users/mccarthy/pkalbers/R-patched/bin/Rscript ../run_psmc.R ../rvage/${CCF}

done
