#!/bin/bash
#$ -N TGP_psm
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e psm_err.txt
#$ -o psm_log.txt
#$ -cwd

cd psm


../rvage/rvage psm -i ../data/1000G.chr${chr}_c.bin -o count.chr${chr} --kList ${frq} -b 5000


