#
# write sharing pair table for Shared Haplotype (approx IBD) detection
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]             # file prefix
shapairs = as.numeric(args[2]) # number of distributed sharing pairs files


load(sprintf("%s.RData", prefix))

G = load.bigmatrix(sprintf("%s.G", prefix))
H = load.bigmatrix(sprintf("%s.H", prefix))
P = load.bigmatrix(sprintf("%s.P", prefix))


dir.create("pairs", showWarnings = F)
setwd("pairs")


cat("Distributing sharing pair tables ... \n")

rng = 1:nrow(FKI)
rng = split(rng, cut(rng, shapairs))

num.pairs = 0

for (k in 1:length(rng)) {
	cat(sprintf("%d of %d  ", k, length(rng)))
	
	pair.name = sprintf("pairs.%s.%04d.RData", prefix, k)
	
	pair = NULL
	
	for (i in rng[[k]]) {
		
		g = parse.sharers(FKI$g.sharer[i])
		h = parse.sharers(FKI$h.sharer[i])
		p = parse.sharers(FKI$p.sharer[i])
		
		g = combn(g, 2, simplify = F)
		h = combn(h, 2, simplify = F)
		p = combn(p, 2, simplify = F)
		
		tmp = data.table(index = FKI$index[i],
										 position = FKI$pos[i],
										 fk = FKI$n.sharer[i],
										 g0 = sapply(g, function(x) x[1]),
										 g1 = sapply(g, function(x) x[2]),
										 h0 = sapply(h, function(x) x[1]),
										 h1 = sapply(h, function(x) x[2]),
										 p0 = sapply(p, function(x) x[1]),
										 p1 = sapply(p, function(x) x[2]))
		
		pair = rbind(pair, tmp)
		
		num.pairs = num.pairs + nrow(tmp)
	}
	
	save(pair, file = pair.name)
	cat("\n")
}

cat("Number of pairs: ", num.pairs, "\n")

