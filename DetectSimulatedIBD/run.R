#
# detect shared haplotype, on genotypes, haplotypes, or phased haplotypes
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]             # file prefix



### detection methods

detect.shap.gen = function(idx, g0, g1, G, POS) {
	g0 = G[, g0]
	g1 = G[, g1]
	
	if (g0[idx] != 1 || g1[idx] != 1) {
		return(c(NA, NA))
	}
	
	brk = abs(g0 - g1)
	brk = which(brk == 2)
	
	lhs = which(brk < idx)
	if (length(lhs) == 0) {
		lhs = 0
	} else {
		lhs = brk[max(lhs)]
	}
	
	rhs = which(brk > idx)
	if (length(rhs) == 0) {
		rhs = length(POS) + 1
	} else {
		rhs = brk[min(rhs)]
	}
	
	return(c(POS[lhs+1], POS[rhs-1]))
}


detect.shap.hap = function(idx, h0, h1, H, POS) {
	h0m = H[, h0]
	h1m = H[, h1]
	
	h0p = if (h0 %% 2 == 0) H[, h0 - 1] else H[, h0 + 1]
	h1p = if (h1 %% 2 == 0) H[, h1 - 1] else H[, h1 + 1]
	
	if (h0m[idx] != h1m[idx] ||
			h0m[idx] == h0p[idx] ||
			h1m[idx] == h1p[idx]) {
		return(c(NA, NA))
	}
	
	brk.m = abs(h0m - h1m)
	brk.p = abs(h0p - h1p)
	
	brk = intersect(which(brk.m == 1), which(brk.p == 1))
	
	lhs = which(brk < idx)
	if (length(lhs) == 0) {
		lhs = 0
	} else {
		lhs = brk[max(lhs)]
	}
	
	rhs = which(brk > idx)
	if (length(rhs) == 0) {
		rhs = length(POS) + 1
	} else {
		rhs = brk[min(rhs)]
	}
	
	return(c(POS[lhs+1], POS[rhs-1]))
}



### execute detections

run.gen = function(pair, G, POS, log) {
	cat("###  Genotypes\n", file = log, append = T)
	
	pair$g.lhs = 0
	pair$g.rhs = 0
	
	for (i in 1:nrow(pair)) {
		if (i %% 1000 == 0) cat(".", file = log, append = T)
		g.lr = detect.shap.gen(pair$index[i], pair$g0[i], pair$g1[i], G, POS)
		pair$g.lhs[i] = g.lr[1]
		pair$g.rhs[i] = g.lr[2]
	}
	cat("\nOK\n", file = log, append = T)
	
	return(pair)
}


run.hap = function(pair, H, POS, log) {
	cat("###  Haplotypes\n", file = log, append = T)
	
	pair$h.lhs = 0
	pair$h.rhs = 0
	
	for (i in 1:nrow(pair)) {
		if (i %% 1000 == 0) cat(".", file = log, append = T)
		h.lr = detect.shap.hap(pair$index[i], pair$h0[i], pair$h1[i], H, POS)
		pair$h.lhs[i] = h.lr[1]
		pair$h.rhs[i] = h.lr[2]
	}
	cat("\nOK\n", file = log, append = T)
	
	return(pair)
}


run.phap = function(pair, P, POS, log) {
	cat("###  Phased haplotypes\n", file = log, append = T)
	
	pair$p.lhs = 0
	pair$p.rhs = 0
	
	for (i in 1:nrow(pair)) {
		if (i %% 1000 == 0) cat(".", file = log, append = T)
		p.lr = detect.shap.hap(pair$index[i], pair$p0[i], pair$p1[i], P, POS)
		pair$p.lhs[i] = p.lr[1]
		pair$p.rhs[i] = p.lr[2]
	}
	cat("\nOK\n", file = log, append = T)
	
	return(pair)
}



### run

load(sprintf("%s.RData", prefix))

G = sprintf("%s.G", prefix)
H = sprintf("%s.H", prefix)
P = sprintf("%s.P", prefix)

G = load.bigmatrix(G)
H = load.bigmatrix(H)
P = load.bigmatrix(P)

G = as.matrix(G)
H = as.matrix(H)
P = as.matrix(P)


setwd("pairs")


pair.files = dir(pattern = sprintf("^pairs\\.%s\\.[0-9]+\\.RData", prefix))


for (pair.file in sample(pair.files)) {
	
	cat("RUN: ", pair.file, "\n")
	
	run = sub("[^0-9]+([0-9]+)[^0-9]+", "\\1", pair.file)
	
	save.file = sprintf("result.%s",  pair.file)
	log.file  = sprintf("log.%s.txt", pair.file)
	
	Sys.sleep(runif(n = 1, min = 1/4, max = 2))
	
	if (file.exists(save.file) || file.exists(log.file)) {
		next
	}
	
	cat("", file = log.file, append = F)
	
	load(pair.file)
	
	if (file.exists(save.file)) next
	
	pair = run.gen(pair, G, POS, log.file)
	
	if (file.exists(save.file)) next
	
	pair = run.hap(pair, H, POS, log.file)
	
	if (file.exists(save.file)) next
	
	pair = run.phap(pair, P, POS, log.file)
	
	if (file.exists(save.file)) next

	cat("Saving results\n", file = log.file, append = T)
	save(pair, file=save.file)
}




