


library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	load(sprintf("pruned.chr%d.pub_pat.RData", i))
	
	
	if (!is.null(pub))
		cat(unique(sample(pub$Position)), sep = "\n", file = sprintf("selected_chr%d.pubmed.txt", i))
	
	if (!is.null(pat))
		cat(unique(sample(pat$Position)), sep = "\n", file = sprintf("selected_chr%d.pathog.txt", i))
	
}









