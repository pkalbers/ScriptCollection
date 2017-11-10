#!/bin/bash

cd psm

for chr in `seq 1 22`
do
	../rvage psm -i ../data/1000G.chr${chr}_c.bin -o count.chr${chr} -k 2 25 -b 5000 
done

