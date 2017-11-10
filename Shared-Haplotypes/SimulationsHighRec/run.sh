#!/bin/bash

# Sequence length:
# 1Mb (1000000)

# Mutation rate per base 
# mutation rate =     1e-08
#   recomb rate =     1e-09
# convert to scrm input:
# mutation rate =     1e-08 * 1000000 * 4 * 10000
#   recomb rate =     1e-09 * 1000000 * 4 * 10000
# expect to see 10 mutations for every recombination event

#./scrm ${NHAP} 1 -t 400 -r 40 1000000 -T -p 15 -SC abs -seed $RANDOM > ${FILE}.out
#./scrm ${NHAP} 1 -t 2000 -r 200 5000000 -T -p 5 -SC abs -seed $RANDOM > ${FILE}.out
#./scrm ${NHAP} 1 -t 8000 -r 800 20000000 -T -p 15 -seed $RANDOM > ${FILE}.out

#grep ";" ${FILE}.out | sed -e 's/\[//g' -e 's/\]/ /g' > ${FILE}.tre &


NHAP=5000
FILE="sim.${NHAP}"


# Simulation

./scrm ${NHAP} 1 -t 4000 -r 4000 10000000 -T -p 10 -seed $RANDOM > ${FILE}.out


# Tree extraction

grep ';' ${FILE}.out | sed -e 's/^/NEWICK TREE: /' > ${FILE}.tre


# IBD detection

cat ${FILE}.tre | ./IBDdetection_naive -l 10000000 -F 1 -m 40000 -t 1 -e 0.00000001 -T 1 -o ${FILE}.tmp
#./IBDdetection_naive -f ${FILE}.tre -l 10000000 -F 1 -m 20000 -t 0 -e 0.00000001 -T 10 -o ${FILE}.tmp


# misc

sed 's/^\([0-9]*\)-\([0-9]*\) \[\([0-9]*\)\,\([0-9]*\)\]$/\1 \2 \3 \4/' ${FILE}.tmp > ${FILE}.ibd &

grep ";" ${FILE}.out | sed -n 's/^\[\([0-9]*\)\].*$/\1/p' > ${FILE}.rec &

grep "segsites:" ${FILE}.out | cut -d' ' -f2- > ${FILE}.seg &

grep "positions:" ${FILE}.out | cut -d' ' -f2- > ${FILE}.pos &

tail -${NHAP} ${FILE}.out > ${FILE}.hap &

wait

rm ${FILE}.tmp


# R parsing

Rscript loadSimData.R $FILE


# clear

#rm ${FILE}.out
#rm ${FILE}.tre
#rm ${FILE}.rec
#rm ${FILE}.seg
#rm ${FILE}.pos
#rm ${FILE}.hap
#rm ${FILE}.tmp
#rm ${FILE}.ibd

