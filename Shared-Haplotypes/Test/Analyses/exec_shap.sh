#!/bin/bash

#
# shared haplotype detectino
#


HAPFILES=( ../GoT2D/converted/GoT2D.*.hap )

for HAPFILE in ${HAPFILES[@]}
do 
	echo $HAPFILE
	
	LEGFILE=${HAPFILE%.*}.legend
	echo $LEGFILE
	
	OUTFILE=${HAPFILE%.*}
	OUTFILE=${OUTFILE%.*}
	OUTFILE=${OUTFILE##*/}
	echo $OUTFILE
	
	./shap \
	-h $HAPFILE \
	-l $LEGFILE \
	-o data.${OUTFILE} \
	-r 0.005 > log.${OUTFILE}.txt 2>&1 &
	
	echo ""
	
done

