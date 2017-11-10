#!/bin/bash
#$ -N vanchkpacks
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e vanchkpacks_err.txt
#$ -o vanchkpacks_log.txt
#$ -cwd

cd packs

while true
do
	sleep $[ ( $RANDOM % 10 ) + 1 ]s

	for file in `ls pack.*.txt | sort -R`
	do
		echo $file
		mv $file _${file}
		/users/mccarthy/pkalbers/miniconda2/bin/python ../get_ibd_result.py _${file} ../vanilla.hdf5 > truth.${file} &
		break
	done
done

