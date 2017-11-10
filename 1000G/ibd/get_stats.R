

pl = NULL
gl = NULL


# FGT

files = dir(pattern = "^result\\.chr([0-9]+)\\.ibd\\.FGT\\.txt$")

for (file in files) {
	print(file)
	CHR = as.numeric(sub("^result\\.chr([0-9]+)\\.ibd\\.FGT\\.txt$", "\\1", file))
	
	d = read.table(file, header = T, stringsAsFactors = F)
	m = read.table(sprintf("1000G.chr%d.marker.txt", CHR), header = T, stringsAsFactors = F)
	
	save(d, m, sprintf("%s.RData", file))
	
	d = d[order(d$Fk), ]
	
	l = d$LHS + 1
	r = d$RHS + 1
	
	PL = m$Position[r] - m$Position[l]
	GL = m$GenDist[r] - m$GenDist[l]
	FK = d$Fk
	
	x = sprintf("%d %d", d$LHS, d$RHS)
	x = which(duplicated(x))
	
	Ntotal = nrow(d)
	Nunique = Ntotal - length(x)
	
	if (length(x) > 0) {
		PL = PL[-x]
		GL = GL[-x]
		FK = FK[-x]
	}
	
	z = sapply(split(PL, FK), fivenum)
	z = t(z)
	z = as.data.frame(z)
	names(z) = c("min","low","med","upp","max")
	z = cbind(chr = CHR, fk = as.numeric(rownames(z)), z)
	
	pl = rbind(pl, z)
	
	
	z = sapply(split(GL, FK), fivenum)
	z = t(z)
	z = as.data.frame(z)
	names(z) = c("min","low","med","upp","max")
	z = cbind(chr = CHR, fk = as.numeric(rownames(z)), z)
	
	gl = rbind(gl, z)
}


save(pl, gl, "stats_fgt.RData")





pl = NULL
gl = NULL


# DGT

files = dir(pattern = "^result\\.chr([0-9]+)\\.ibd\\.DGT\\.txt$")

for (file in files) {
	print(file)
	CHR = as.numeric(sub("^result\\.chr([0-9]+)\\.ibd\\.FGT\\.txt$", "\\1", file))
	
	d = read.table(file, header = T, stringsAsFactors = F)
	m = read.table(sprintf("1000G.chr%d.marker.txt", CHR), header = T, stringsAsFactors = F)
	
	save(d, m, sprintf("%s.RData", file))
	
	d = d[order(d$Fk), ]
	
	l = d$LHS + 1
	r = d$RHS + 1
	
	PL = m$Position[r] - m$Position[l]
	GL = m$GenDist[r] - m$GenDist[l]
	FK = d$Fk
	
	x = sprintf("%d %d", d$LHS, d$RHS)
	x = which(duplicated(x))
	
	Ntotal = nrow(d)
	Nunique = Ntotal - length(x)
	
	if (length(x) > 0) {
		PL = PL[-x]
		GL = GL[-x]
		FK = FK[-x]
	}
	
	z = sapply(split(PL, FK), fivenum)
	z = t(z)
	z = as.data.frame(z)
	names(z) = c("min","low","med","upp","max")
	z = cbind(chr = CHR, fk = as.numeric(rownames(z)), z)
	
	pl = rbind(pl, z)
	
	
	z = sapply(split(GL, FK), fivenum)
	z = t(z)
	z = as.data.frame(z)
	names(z) = c("min","low","med","upp","max")
	z = cbind(chr = CHR, fk = as.numeric(rownames(z)), z)
	
	gl = rbind(gl, z)
}


save(pl, gl, "stats_dgt.RData")







