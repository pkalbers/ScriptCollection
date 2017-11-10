#!/bin/bash

# Sequence length:
# 50Mb (50000000)

# mutation rate =     1e-08
#   recomb rate =     1e-08
# convert to scrm input:
# mutation rate =     1e-08 * 50000000 * 4 * 10000
#   recomb rate =     1e-08 * 50000000 * 4 * 10000
# expect to see 1 mutations for every recombination event



NHAP=2000
FILE="sim.${NHAP}"


# Simulation

./scrm ${NHAP} 1 -t 20000 -r 20000 50000000 -T -p 10 -seed $RANDOM > ${FILE}.out


grep "segsites:" ${FILE}.out | cut -d' ' -f2- > ${FILE}.seg &

grep "positions:" ${FILE}.out | cut -d' ' -f2- > ${FILE}.pos &

tail -${NHAP} ${FILE}.out > ${FILE}.hap &

wait


