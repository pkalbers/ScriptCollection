#!/bin/bash

CORES="12"

GENFILE=`ls *.gen`
SAMFILE=`ls *.sample`
PHSFILE=`basename $GENFILE`
PHSFILE=${GENFILE%.gen}
PHSFILE="phased.${PHSFILE}"

echo "GEN file:      ${GENFILE}"
echo "Sample file:   ${SAMFILE}"
echo "Output prefix: ${PHSFILE}"

../shapeit -G $GENFILE $SAMFILE -O $PHSFILE -T $CORES

