#!/bin/bash

#
# recombination rate:  1e-8 * 4 * 10000 = 0.0004
# mutation rate:     1.2e-8 * 4 * 10000 = 0.00048
# chromsome length:  100Mb = 100000000
#

PREFIX="history"  # output filename prefix
NUMSAMPLE="5000"   # sample size
NUMMARKER="100000000"   # number of loci
RECRATE="0.0004"   # recombination rate, scaled
MUTRATE="0.00048"   # mutation rate, scaled

FKMAX="25"
NUMPAIRFILES="64"

python3 simulate.py -o ${PREFIX} -s ${NUMSAMPLE} -n ${NUMMARKER} -r ${RECRATE} -m ${MUTRATE}

python3 get_positions.py -i ${PREFIX}

python3 get_haplotypes.py -i ${PREFIX} -t False

Rscript load.R ${FKMAX} ${PREFIX} ${NUMPAIRFILES}
