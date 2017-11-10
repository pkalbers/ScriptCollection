


library(data.table)


files = dir(pattern = ".+marker.txt$")

M = list()

for (file in files) {
	
	chr = as.numeric(sub(".+chr([0-9]+)_c\\..+", "\\1", file))
	
	cat(chr)
	
	M[[chr]] = fread(file, header = T, stringsAsFactors = F)
	
}


files = dir(pattern = "random.+")

res = NULL

for (file in files) {
	
	chr = as.numeric(sub(".+chr([0-9]+)\\..+", "\\1", file))
	cat(".")
	
	tmp = read.table(file, header = T, stringsAsFactors = F)
	
	if (nrow(tmp) == 0) next
	
	tmp = as.data.table(tmp)
	tmp = cbind(tmp, chr = chr)
	
	res = rbind(res, tmp)
	
}


files = dir(pattern = "target.+")

tar = NULL

for (file in files) {
	
	chr = as.numeric(sub(".+chr([0-9]+)\\..+", "\\1", file))
	cat(".")
	
	tmp = read.table(file, header = T, stringsAsFactors = F)
	
	if (nrow(tmp) == 0) next
	
	tmp = as.data.table(tmp)
	tmp = cbind(tmp, chr = chr)
	
	tar = rbind(tar, tmp)
	
}




d = rbind(res, tar)

x = data.frame(chr = d$chr,  mid = d$MarkerID)

out = NULL

for (i in 1:nrow(x)) {
	cat(".")
	z = which(M[[ x$chr[i] ]]$MarkerID == x$mid[i])
	
	out = rbind(out, M[[ x$chr[i] ]][z, ])
}


out$A0 = sub("^([A-Z]),([A-Z])$", "\\1", out$Alleles)
out$A1 = sub("^([A-Z]),([A-Z])$", "\\2", out$Alleles)


rs = data.table(CHR = out$Chromosome, POS = out$Position, A0 = out$A0, A1 = out$A1, RSID = out$Label)

rs$RSID = sub(".+_(rs[0-9]+)$", "\\1", rs$RSID)

rs = rs[grep("^rs.+", rs$RSID),]

rs = rs[order(rs$CHR, rs$POS), ]

vars = rs


save(vars, file = "UKBB_age_vars.RData")

write.table(vars, file = "UKBB_age_vars.txt", quote = F, row.names = F, col.names = T)








