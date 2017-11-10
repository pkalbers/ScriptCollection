#!/bin/bash

CORES=${1:-4}

GENFILE=`ls *.gen`
SAMFILE=`ls *.sample`
PHSFILE=`basename $GENFILE`
PHSFILE=${GENFILE%.gen}
PHSFILE="phased.${PHSFILE}"

echo ""
echo "GEN file:    ${GENFILE}"
echo "SAMPLE file: ${SAMFILE}"
echo "# cores:     ${CORES}"
echo ""

./shapeit -G ${GENFILE} ${SAMFILE} -O ${PHSFILE} -T ${CORES}


