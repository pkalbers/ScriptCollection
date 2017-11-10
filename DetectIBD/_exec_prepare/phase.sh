#!/bin/bash

CORES="12"

GENFILE=`ls *.gen`
SAMFILE=`ls *.sample`
PHSFILE=`basename $GENFILE`
PHSFILE=${GENFILE%.gen}
PHSFILE="phased.${PHSFILE}"

MAPFILE="../shapeit.3col.genetic_map_GRCh37_chr20.txt"

echo "GEN file:      ${GENFILE}"
echo "Sample file:   ${SAMFILE}"
echo "Output prefix: ${MAPFILE}"
echo "Output prefix: ${PHSFILE}"

../shapeit -G $GENFILE $SAMFILE -O $PHSFILE -M $MAPFILE -T $CORES

