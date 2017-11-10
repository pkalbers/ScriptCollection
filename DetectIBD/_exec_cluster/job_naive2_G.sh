#!/bin/bash
#$ -N shap_G_detect
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _shap_G_detect_err.txt
#$ -o _shap_G_detect_log.txt
#$ -cwd

~/R-patched/bin/Rscript /well/pkalbers/detect_IBD/run_naive.R OutOfAfricaHapMap20.GenErr_1000G G
