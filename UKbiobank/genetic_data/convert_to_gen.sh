#!/bin/bash


for CHR in `seq 1 22`
do
	mkdir chr${CHR}
done

for CHR in `seq 1 22`
do
	BIN="chr${CHR}impv1.bgen"
	GEN="chr${CHR}impv1.gen"
	DIR="chr${CHR}"
	OUT="qctool.chr${CHR}impv1.txt"
	LOG="qctool.chr${CHR}impv1.log"

	echo
	echo $BIN
	echo $GEN
	echo $DIR
	echo $OUT
	echo $LOG
	echo

	nohup ./qctool -g ${BIN} -og ./${DIR}/${GEN} -log ./${DIR}/${LOG} > ./${DIR}/${OUT} &

done

