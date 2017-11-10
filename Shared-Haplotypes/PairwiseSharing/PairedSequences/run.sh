#!/bin/bash

./ship -i ../merged_related.vcf.gz -f 5 -s ../merged_related.panel -m ../genetic_map_GRCh37_chr20.txt --threads 6 --remove_unknown_markers --remove_multiallelic -sub ../Trios/PUR.sub_panel -o sequence_PUR

./ship -i ../merged_related.vcf.gz -f 5 -s ../merged_related.panel -m ../genetic_map_GRCh37_chr20.txt --threads 6 --remove_unknown_markers --remove_multiallelic -sub ../Trios/YRI.sub_panel -o sequence_YRI

./ship -i ../merged_related.vcf.gz -f 5 -s ../merged_related.panel -m ../genetic_map_GRCh37_chr20.txt --threads 6 --remove_unknown_markers --remove_multiallelic -sub ../Trios/KHV.sub_panel -o sequence_KHV

./ship -i ../merged_related.vcf.gz -f 5 -s ../merged_related.panel -m ../genetic_map_GRCh37_chr20.txt --threads 6 --remove_unknown_markers --remove_multiallelic -sub ../Trios/MXL.sub_panel -o sequence_MXL
