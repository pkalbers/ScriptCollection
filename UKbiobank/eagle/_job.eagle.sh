#!/bin/bash
#$ -N EagleUKB
#$ -P mccarthy.prjc
#$ -q himem.qh
#$ -e _err.eagle.txt
#$ -o _log.eagle.txt
#$ -cwd
#$ -pe shmem 48


./eagle \
--bfile chr7impv1.plink \
--geneticMapFile eagle_genetic_map_hg19_withX.txt.gz \
--outPrefix chr7impv1.eagle \
--numThreads 48 \
--chrom 7

