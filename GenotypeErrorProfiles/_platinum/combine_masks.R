
library(data.table)


A = "./20140520.strict_mask.autosomes.bed"
B = "./ConfidentRegions.bed"


a = fread(A, header = F, stringsAsFactors = F)
b = fread(B, header = F, stringsAsFactors = F)


a = split(a, a$V1)
b = split(b, b$V1)


tag = intersect(names(a), names(b))

out = NULL

for (chr in tag) {
	print(chr)
	
	ma = a[[chr]]
	mb = b[[chr]]
	
	end = max(c(ma$V3, mb$V3))
	
	xa = rep(F, end)
	xb = rep(F, end)
	
	for (i in 1:nrow(ma)) {
		rng = (ma$V2[i]):(ma$V3[i] - 1)
		xa[rng] = T
	}
	
	for (i in 1:nrow(mb)) {
		rng = (mb$V2[i]):(mb$V3[i] - 1)
		xb[rng] = T
	}
	
	x = xa & xb
	
	z = which(x[1:(end-1)] != x[2:(end)])
	z = z + 1
	
	out = rbind(out, data.table(chr = chr, 
															beg = z[(1:length(z)) %% 2 != 0],
															end = z[(1:length(z)) %% 2 == 0]))
	
}



write.table(out, file = "CombinedRegions.bed", quote = F, row.names = F, col.names = F)


