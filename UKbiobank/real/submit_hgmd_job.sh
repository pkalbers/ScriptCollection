#!/bin/bash


for CHR in `seq 1 22`
do
	FILE="./hgmd/hgmd_chr${CHR}.txt"

	while read LINE
	do
		REF=`echo $LINE | cut -d'|' -f1`
		POS=`echo $LINE | cut -d'|' -f2`

		echo "chr ${CHR} ref ${REF} pos ${POS}"

		qsub -v ref=${REF},chr=${CHR},pos=${POS} ./hgmd_job.sh

		echo

	done < $FILE

done

