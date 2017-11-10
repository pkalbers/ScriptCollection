#!/bin/bash
#$ -N ESH
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e ./logs/err.ESH.txt
#$ -o ./logs/err.ESH.txt
#$ -cwd


export PATH=/apps/well/R/3.1.3/bin:$PATH

Rscript ./runESH.R data.sim.2000.RData result_ESH > ./logs/log.ESH.${RANDOM}.txt
