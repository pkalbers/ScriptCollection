#!/bin/bash

VCF=$1
TAG=$2
OUT="subset.${VCF}.vcf"

> $OUT

gunzip -c $VCF | while read line
do

	PRE=`echo "${line}" | cut -c 1`

	if [ "${PRE}" = "#" ]
	then
		echo "HEAD"
		#echo $line >> $OUT
		continue
	fi

	SUB=`echo "${line}" | cut -f 2 -d$'\t' `

	echo $SUB

	while read tag
	do
		if [ $SUB -eq $tag ]
		then
			echo "^^^^^^^^^"
			echo $line >> $OUT
			break
		fi

	done < $TAG


done

