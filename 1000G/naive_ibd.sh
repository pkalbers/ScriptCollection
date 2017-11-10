#!/bin/bash

cd length

for chr in `seq 1 22`
do
	../rvage ibd -i ../data/1000G.chr${chr}_c.bin -o length2.chr${chr} -k 2 26 --randomPairs 100000 -b 5000 -m fgt -t 4
	../rvage ibd -i ../data/1000G.chr${chr}_c.bin -o length2.chr${chr} -k 2 26 --randomPairs 100000 -b 5000 -m dhg -t 4
done

