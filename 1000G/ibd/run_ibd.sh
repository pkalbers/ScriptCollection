#!/bin/bash  

for chr in `seq 1 22`
do

../../rvage/rvage ibd -i ../data/1000G.chr${chr}.bin -o result.chr${chr} -k 2 25 -b 5000 -m fgt -t 60

done

