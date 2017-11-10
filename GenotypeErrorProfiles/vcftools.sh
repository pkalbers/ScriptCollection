#!/bin/bash

FILE=$1
INDV=$2

echo
echo "Input file:    ${FILE}"
echo

vcftools --gzvcf $FILE \
	--indv ${INDV} --remove-indels --recode \
	--out ${INDV}_${FILE}

echo
echo "DONE"
echo
