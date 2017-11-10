#!/bin/bash

HFILE=$1
CORES=$2

echo ""
echo "History file:    ${HFILE}"
echo "Number of cores: ${CORES}"
echo ""

for PFILE in truth.*.txt
do
	HBASE=`basename ${HFILE}`
	PBASE=`basename ${PFILE}`

	THFILE="/tmp/IBD_DETECT_IN_PAIRS.${PBASE}.${HBASE}"
	cp $HFILE $THFILE

	LFILE="log.age.${PBASE}"

	echo ""
	echo "History file:  ${THFILE}"
	echo "Sharing Pairs: ${PFILE}"
	echo "Log file:      ${LFILE}"
	echo ""

	python ../detect_age.py -i ${THFILE%.*} -p $PFILE > $LFILE &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $CORES )) ### max number of parallel processes
	do
        sleep 5
        JOBLIST=($(jobs -p))
    done
done

wait

rm /tmp/IBD_DETECT_IN_PAIRS.*
