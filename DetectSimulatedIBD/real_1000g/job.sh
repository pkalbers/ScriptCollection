#!/bin/bash
#$ -N shap_detect
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _shap_detect_err.txt
#$ -o _shap_detect_log.txt
#$ -cwd
#$ -pe shmem 4

~/R-patched/bin/Rscript run.R
