#!/bin/bash


for file in `ls pack.*`
do 
	echo $file
	
	python get_ibd_result.py ${file} vanilla.hdf5 > truth.${file} &
	
done
wait


