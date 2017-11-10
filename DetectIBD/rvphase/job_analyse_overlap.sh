#!/bin/bash
#$ -N rvphase_over
#$ -P mccarthy.prjc
#$ -q long.qc
#$ -e _rvphase_overlap_err.txt
#$ -o _rvphase_overlap_log.txt
#$ -cwd

export PATH=/opt/well/R/R-3.0.2/bin:$PATH
export PATH=/apps/well/R/3.1.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/well/R/R-3.0.2/lib64/R/lib:/opt/well/libpng/1.5.18/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/apps/well/R/3.1.0/lib64/R/lib:$LD_LIBRARY_PATH

cd ./overlap/

Rscript ../analyse_overlap.R

