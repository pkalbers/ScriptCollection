

library(ggplot2)
library(grid)
library(bigmemory)
library(data.table)


data.file = "data.sim.5000.RData"

load(data.file)


G = load.bigmatrix(gen.file)

l = nrow(M)

line = vector(mode = "character", length = l)

for (i in 1:l) {
	cat(i, "of", l, "\n")
	
	chr = 1
	snp = sprintf("SNP_%d", i)
	pos = M$position[i]
	al0 = sample(c("A", "C", "G", "T"), 1)
	al1 = sample(c("A", "C", "G", "T"), 1) 
	while (al1 == al0) {
		al1 = sample(c("A", "C", "G", "T"), 1) 
	}
	
	gen = sapply(G[i, ], function(g) {
		if (g == 0) return("1 0 0")
		if (g == 1) return("0 1 0")
		if (g == 2) return("0 0 1")
	})
	
	line[i] = paste(chr, snp, pos, al0, al1, paste(gen, collapse = " "))
}


writeLines(line, "./phasing_compare/sim.5000.gen")



sam = data.table(id1 = c("ID_1", "0", as.vector(J)), id1 = c("ID_2", "0", as.vector(J)), miss = c("missing", "0", rep("0", length(J))))

write.table(sam, file = "./phasing_compare/sim.5000.sample", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)












