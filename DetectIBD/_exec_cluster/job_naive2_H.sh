#!/bin/bash
#$ -N shap_H_detect
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _shap_H_detect_err.txt
#$ -o _shap_H_detect_log.txt
#$ -cwd
#$ -pe shmem 2

~/R-patched/bin/Rscript /well/pkalbers/detect_IBD/run_naive.R OutOfAfricaHapMap20.GenErr_1000G H
