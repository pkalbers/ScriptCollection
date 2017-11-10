

library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	tmp = read.table(sprintf("chr%d_imp_hg38.vcf.gz", i), header = F, stringsAsFactors = F)
	
	vcf = tstrsplit(tmp[[3]], split=":", fixed=T, type.convert = T)
	vcf = as.data.table(vcf)
	names(vcf) = c("chr", "pos", "a0", "a1")
	
	vcf$POS38 = tmp[[2]]
	
	
	vep = read.table(sprintf("chr%d_imp_hg38_vep_anns.txt.gz", i), header = F, stringsAsFactors = F)
	#ext = tmp[, 14]
	#tmp = tmp[, 1:13]
	names(vep) = c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type",
								 "Consequence","cDNA_position","CDS_position","Protein_position",
								 "Amino_acids","Codons","Existing_variation","Extras")
	vep = as.data.table(vep)
	
	
	save(vcf, vep, file = sprintf("chr%d.vcf_vep.RData", i))
}


