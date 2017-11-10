#!/bin/bash
#$ -N rvphase_shap
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e _rvphase_shap_err.txt
#$ -o _rvphase_shap_log.txt
#$ -cwd
#$ -pe shmem 2

export PATH=/opt/well/R/R-3.0.2/bin:$PATH
export PATH=/apps/well/R/3.1.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/well/R/R-3.0.2/lib64/R/lib:/opt/well/libpng/1.5.18/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/apps/well/R/3.1.0/lib64/R/lib:$LD_LIBRARY_PATH

cd ./shap/

Rscript ../analyse_shap.R

