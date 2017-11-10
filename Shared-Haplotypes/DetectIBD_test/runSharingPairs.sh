#!/bin/bash

# detect history file
HFILE=`ls *.hdf5`

echo "Using ${HFILE}"

for PFILE in sharingpairs.*
do
	THFILE=`basename ${PFILE}`
	THFILE="/tmp/${THFILE}.${HFILE}"
	cp $HFILE $THFILE

	LFILE="log.${PFILE}"

	echo "History file:  ${HFILE}"
	echo "Sharing Pairs: ${PFILE}"
	echo "Log file:      ${LFILE}"

	python get_ibd_segments.py $THFILE $PFILE > $LFILE &

	JOBLIST=($(jobs -p))
	while (( ${#JOBLIST[*]} >= 20 )) ### max number of parallel processes
	do
        sleep 30
        JOBLIST=($(jobs -p))
    done
done

wait
