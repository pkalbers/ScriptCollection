#!/bin/bash

PREFIX=$1
CORES=$2
NJOBS=$3
SUFFIX="RData"

ROOT=`pwd`

for FILE in `find -iname "${PREFIX}*${SUFFIX}"`;
do

	cd `dirname $FILE`

	BASE=`basename $FILE`
	BASE=${BASE%.RData}

	echo "###"
	echo $BASE
	echo "###"


	GENFILE=`ls *.gen`
	SAMFILE=`ls *.sample`
	PHSFILE=`basename $GENFILE`
	PHSFILE=${GENFILE%.gen}
	PHSFILE="phased.${PHSFILE}"

	echo "GEN file:      ${GENFILE}"
	echo "Sample file:   ${SAMFILE}"
	echo "Output prefix: ${PHSFILE}"

	../shapeit -G $GENFILE $SAMFILE -O $PHSFILE -T $CORES &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $NJOBS ))
	do
        sleep 30
        JOBLIST=($(jobs -p))
    done

	echo ""

	cd $ROOT

done

wait
