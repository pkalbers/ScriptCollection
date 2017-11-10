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

	LFILE="log.${PBASE}"

	echo ""
	echo "History file:  ${THFILE}"
	echo "Sharing Pairs: ${PFILE}"
	echo "Log file:      ${LFILE}"
	echo ""

	python ../detect_ibd.py -i ${THFILE%.*} -p $PFILE > $LFILE &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $CORES )) ### max number of parallel processes
	do
        sleep 30
        JOBLIST=($(jobs -p))
    done
done

wait

rm /tmp/IBD_DETECT_IN_PAIRS.*
