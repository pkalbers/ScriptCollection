#!/bin/bash

PRX="truth"
OUT="simres.dis1000.txt"
FST="${PRX}.pack.0001.*"

head -n 1 `ls $FST` > $OUT

for file in `ls ${PRX}*`
do 
	echo ${file}
	#head -n 200001 $file | tail -n 200000 >> $OUT
	tail -n +2 ${file} >> $OUT
done


