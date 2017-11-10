#!/bin/bash


for file in `ls pack.*`
do 
	echo $file
	
	/flash/compK000/pkalbers/miniconda2/bin/python get_ibd_result.py ${file} vanilla.hdf5 > truth.${file} &
	
done
wait


