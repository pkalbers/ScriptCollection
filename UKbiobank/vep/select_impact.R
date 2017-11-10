


library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	load(sprintf("pruned.chr%d.vep.RData", i))
	
	tmp = split(vep, vep$impact)
	tmp = lapply(tmp, function(x) { if (nrow(x) > 5000) x[sample(1:nrow(x), size = min(5000, nrow(x)), replace = F), ] else x })
	
	for (tag in names(tmp)) {
		
		sub = tmp[[tag]]
		
		if (nrow(sub) < 10) {
			stop("???")
		}
		
		cat(sample(sub$Position), sep = "\n", file = sprintf("selected_chr%d.impact.%s.txt", i, tag))
		
	}
	
}




