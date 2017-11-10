#!/bin/bash
#$ -N shap_hmm_detect_err
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _shap_hmm_detect_err.txt
#$ -o _shap_hmm_detect_log.txt
#$ -cwd

~/R-patched/bin/Rscript /well/pkalbers/detect_IBD/run_hmm.R OutOfAfricaHapMap20.GenErr_1000G
