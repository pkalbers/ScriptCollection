#!/bin/bash

HFILE=$1
CORES=$2

echo ""
echo "History file:    ${HFILE}"
echo "Number of cores: ${CORES}"
echo ""

for PFILE in pairs.*.txt
do
	HBASE=`basename ${HFILE}`
	PBASE=`basename ${PFILE}`
	LFILE="log.${PBASE}"

	if [ -e $LFILE ]
	then
		echo "  Already processed:  ${PFILE}"
		echo ""
		continue
	fi

	THFILE="/tmp/IBD_DETECT_LOCAL.${PBASE}.${HBASE}"
	cp $HFILE $THFILE

	echo "History file:  ${THFILE}"
	echo "Sharing Pairs: ${PFILE}"
	echo "Log file:      ${LFILE}"
	echo ""

	python ../_exec_msprime/msprime.local_ibd.py -i $THFILE -p $PFILE > $LFILE &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= $CORES )) ### max number of parallel processes
	do
        sleep 30
        JOBLIST=($(jobs -p))
    done
done

wait
