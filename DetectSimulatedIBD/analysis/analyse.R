#
# analyse all results
#


library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


get.fki = function(file) { load(file); FKI }
FKI = list()
FKI[["1000G"]] = get.fki("../generror_1000g/history.generror_1000g.RData")
FKI[["AFFYMETRIXAXIOM"]] = get.fki("../generror_affymetrixaxiom/history.generror_affymetrixaxiom.RData")
FKI[["SIMULATED"]] = get.fki("../history.RData")


load("./result.pair.match.RData")

D = split(match, match$setup)


G = list()
H = list()
P = list()

setwd("../generror_1000g/")
G[["1000G"]] = attach.big.matrix(dget("_bigmatrix.history.generror_1000g.G.desc"))
H[["1000G"]] = attach.big.matrix(dget("_bigmatrix.history.generror_1000g.H.desc"))
P[["1000G"]] = attach.big.matrix(dget("_bigmatrix.history.generror_1000g.P.desc"))

setwd("../generror_affymetrixaxiom/")
G[["AFFYMETRIXAXIOM"]] = attach.big.matrix(dget("_bigmatrix.history.generror_affymetrixaxiom.G.desc"))
H[["AFFYMETRIXAXIOM"]] = attach.big.matrix(dget("_bigmatrix.history.generror_affymetrixaxiom.H.desc"))
P[["AFFYMETRIXAXIOM"]] = attach.big.matrix(dget("_bigmatrix.history.generror_affymetrixaxiom.P.desc"))

setwd("../")
G[["SIMULATED"]] = attach.big.matrix(dget("_bigmatrix.history.G.desc"))
H[["SIMULATED"]] = attach.big.matrix(dget("_bigmatrix.history.H.desc"))
P[["SIMULATED"]] = attach.big.matrix(dget("_bigmatrix.history.P.desc"))

setwd("./analysis/")


AC = list()
AF = list()
for (tag in names(H)) {
	AC[[tag]] = rowSums(as.matrix(H[[tag]]))
	AF[[tag]] = AC[[tag]] / ncol(H[[tag]])
}




d = rbind(data.table(setup = equal.match$setup, side = "LHS", fk = equal.match$fk, brk = equal.match$g.brk.lhs, type = "Genotype",         dist = equal.match$g.brk.dist.lhs / equal.match$t.brk.dist.lhs * -1),
					data.table(setup = equal.match$setup, side = "RHS", fk = equal.match$fk, brk = equal.match$g.brk.rhs, type = "Genotype",         dist = equal.match$g.brk.dist.rhs / equal.match$t.brk.dist.rhs ),
					data.table(setup = equal.match$setup, side = "LHS", fk = equal.match$fk, brk = equal.match$h.brk.lhs, type = "Haplotype",        dist = equal.match$h.brk.dist.lhs / equal.match$t.brk.dist.lhs * -1),
					data.table(setup = equal.match$setup, side = "RHS", fk = equal.match$fk, brk = equal.match$h.brk.rhs, type = "Haplotype",        dist = equal.match$h.brk.dist.rhs / equal.match$t.brk.dist.rhs ),
					data.table(setup = equal.match$setup, side = "LHS", fk = equal.match$fk, brk = equal.match$p.brk.lhs, type = "Phased haplotype", dist = equal.match$h.brk.dist.lhs / equal.match$t.brk.dist.lhs * -1),
					data.table(setup = equal.match$setup, side = "RHS", fk = equal.match$fk, brk = equal.match$p.brk.rhs, type = "Phased haplotype", dist = equal.match$h.brk.dist.rhs / equal.match$t.brk.dist.rhs ))




for (tag in names(G)) {
	
	sub = which(d$setup == tag)
	sub = d[sub, ]
	
	g = which(sub$type == "Genotype")
	g = sub[g, ]
	
	break
	
	gu = which(abs(g$dist) < 1)
	gu = g[gu, ]
	
	brk = unique(gu$brk)
	
	G[[tag]][, ]
}


AC = rowSums(H)


z = which((match$g.brk.dist.lhs / match$t.brk.dist.lhs) < 1)






d = rbind(data.table(setup = match$setup, side = "LHS", fk = match$fk, brk = match$g.brk.lhs, type = "Genotype",         dist = match$g.brk.dist.lhs / match$t.brk.dist.lhs * -1),
					data.table(setup = match$setup, side = "RHS", fk = match$fk, brk = match$g.brk.rhs, type = "Genotype",         dist = match$g.brk.dist.rhs / match$t.brk.dist.rhs ),
					data.table(setup = match$setup, side = "LHS", fk = match$fk, brk = match$h.brk.lhs, type = "Haplotype",        dist = match$h.brk.dist.lhs / match$t.brk.dist.lhs * -1),
					data.table(setup = match$setup, side = "RHS", fk = match$fk, brk = match$h.brk.rhs, type = "Haplotype",        dist = match$h.brk.dist.rhs / match$t.brk.dist.rhs ),
					data.table(setup = match$setup, side = "LHS", fk = match$fk, brk = match$p.brk.lhs, type = "Phased haplotype", dist = match$h.brk.dist.lhs / match$t.brk.dist.lhs * -1),
					data.table(setup = match$setup, side = "RHS", fk = match$fk, brk = match$p.brk.rhs, type = "Phased haplotype", dist = match$h.brk.dist.rhs / match$t.brk.dist.rhs ))

d = split(d, list(d$type, d$setup))

z = which((d$g.brk.dist.lhs / match$t.brk.dist.lhs) < 1)



ggplot(data = d[sample(nrow(d), 1e5), ]) +
	facet_grid(setup~side) +
	geom_freqpoly(aes(dist, ..density.., colour = type), binwidth = 10)


z = which(abs(d$dist) < 1)
p = d[z,]
p = split(p, list(p$setup, p$type))
p = lapply(p, function(x) {
	nf = which(!is.finite(x$dist))
	if (length(nf) > 0) x$dist[nf] = 0
	lhs = density(x$dist, bw = 0.005, n = 1024*10, from=-2, to=0)
	lhs = data.table(setup=x$setup[1], type=x$type[1], x = lhs$x, y = lhs$y / sum(lhs$y))
	rhs = density(x$dist, bw = 0.005, n = 1024*10, from=0, to=2)
	rhs = data.table(setup=x$setup[1], type=x$type[1], x = rhs$x, y = rhs$y / sum(rhs$y))
	rbind(lhs, rhs)
})
p = rbindlist(p)


ggplot(data=p) + 
	facet_wrap(~setup, ncol=1, scales = "free_y") + 
	geom_vline(xintercept=c(-1, 0, 1)) +
	geom_line(aes(x=x, y=y, colour=type), alpha=0.75, size=3/4) + 
	scale_x_continuous(breaks = -2:2) +
	theme_bw() +
	xlab("Relative physical distance") +
	ylab("Density")








d = rbind(data.table(setup = equal.match$setup, side="LHS", type = "Genotype",         det = log10(-1 * equal.match$g.brk.dist.lhs + 10), tru = log10(-1 * equal.match$t.brk.dist.lhs + 10)),
					data.table(setup = equal.match$setup, side="RHS", type = "Genotype",         det = log10(     equal.match$g.brk.dist.rhs + 10), tru = log10(     equal.match$t.brk.dist.rhs + 10)),
					data.table(setup = equal.match$setup, side="LHS", type = "Haplotype",        det = log10(-1 * equal.match$h.brk.dist.lhs + 10), tru = log10(-1 * equal.match$t.brk.dist.lhs + 10)),
					data.table(setup = equal.match$setup, side="RHS", type = "Haplotype",        det = log10(     equal.match$h.brk.dist.rhs + 10), tru = log10(     equal.match$t.brk.dist.rhs + 10)),
					data.table(setup = equal.match$setup, side="LHS", type = "Phased haplotype", det = log10(-1 * equal.match$p.brk.dist.lhs + 10), tru = log10(-1 * equal.match$t.brk.dist.lhs + 10)),
					data.table(setup = equal.match$setup, side="RHS", type = "Phased haplotype", det = log10(     equal.match$p.brk.dist.rhs + 10), tru = log10(     equal.match$t.brk.dist.rhs + 10)))


d = lapply(split(d, list(d$setup, d$type, d$side)), function(x) {
	brk = seq(log10(0 + 10), log10(100e6 + 10), length.out = 201)
	det = as.numeric(cut(x$det, breaks = brk, include.lowest = T))
	tru = as.numeric(cut(x$tru, breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = 200, ncol = 200)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	#mat = mat / max(mat)
	mat = melt(mat, value.name = "x")
	mat$setup = x$setup[1]
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	as.data.table(mat)
})
p = rbindlist(d)


p$Var1 = p$Var1 + 1
x = (p$side == "LHS")
p$Var1[ x] = (p$Var1[ x]) * -1


xtik = c(-202, -150, -100, -50, 50, 100, 150, 202); xtag =  c("100Mb", "1Mb", "10Kb", "100b", "100b", "10Kb", "1Mb", "100Mb")
ytik = c(1, 50, 100, 150, 200); ytag = c("1b", "100b", "10Kb", "1Mb", "100Mb")

#ggplot(data = p[which(p$setup=="SIMULATED" & p$type == "Phased haplotype"), ]) + 
gg = ggplot(data = p) + 
	facet_grid(setup~type) + 
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(-1, 1), intercept = c(-2, -2), colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6) +
	#scale_fill_gradient(low = "darkblue", high = "firebrick1", na.value = "grey85") +
	#geom_vline(xintercept = 0, colour="white") +
	annotate("text", label = "LHS", x = -20, y = 190, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x =  20, y = 190, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = xtik, labels = xtag) +
	scale_y_continuous(expand = c(0, 0), breaks = ytik, labels = ytag) +
	theme_classic() +
	theme(aspect.ratio = 1/2) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance (IBD)")

ggsave(filename = "_plot.breaks-truth-inferred.pdf", plot = gg, width = 15, height = 15)









d = rbind(data.table(setup = data$setup, type = "Genotype",         fk = data$fk, Length = data$g.size),
					data.table(setup = data$setup, type = "Haplotype",        fk = data$fk, Length = data$h.size),
					data.table(setup = data$setup, type = "Phased haplotype", fk = data$fk, Length = data$p.size))

d$fk = as.factor(d$fk)
d$type = as.factor(d$type)


ggplot(data = d[sample(nrow(d), 5000000), ]) +
	facet_grid(setup~fk) + 
	geom_boxplot(aes(x = type, y = Length), outlier.shape = NA)










#
# re-check genotyping error
# after assigned true genotypes (TG) and error genotypes (EG)
#
m = matrix(0, 3, 3, dimnames = list(as.character(0:2), as.character(0:2)))

for (i in 1:ncol(TG)) {
	cat(i, "of", ncol(TG), "\n") 
	tg = TG[, i]
	eg = EG[, i]
	r0 = which(tg == 0)
	r1 = which(tg == 1)
	r2 = which(tg == 2)
	
	r0c0 = which(eg[r0] == 0)
	r0c1 = which(eg[r0] == 1)
	r0c2 = which(eg[r0] == 2)
	
	r1c0 = which(eg[r1] == 0)
	r1c1 = which(eg[r1] == 1)
	r1c2 = which(eg[r1] == 2)
	
	r2c0 = which(eg[r2] == 0)
	r2c1 = which(eg[r2] == 1)
	r2c2 = which(eg[r2] == 2)
	
	m["0", "0"] = m["0", "0"] + length(r0c0)
	m["0", "1"] = m["0", "1"] + length(r0c1)
	m["0", "2"] = m["0", "2"] + length(r0c2)
	
	m["1", "0"] = m["1", "0"] + length(r1c0)
	m["1", "1"] = m["1", "1"] + length(r1c1)
	m["1", "2"] = m["1", "2"] + length(r1c2)
	
	m["2", "0"] = m["2", "0"] + length(r2c0)
	m["2", "1"] = m["2", "1"] + length(r2c1)
	m["2", "2"] = m["2", "2"] + length(r2c2)
}





#
# explore simulared error
#

setwd("../")

H = load.bigmatrix(sprintf("history.H"))
G = load.bigmatrix(sprintf("history.G"))
P = load.bigmatrix(sprintf("history.P"))

sim = which(match$setup == "SIMULATED")
sim = match[sim, ]

ge = which(sim$t.brk.lhs < sim$g.brk.lhs | sim$t.brk.rhs > sim$g.brk.rhs)
ge = sim[ge, ]

he = which(sim$t.brk.lhs < sim$h.brk.lhs | sim$t.brk.rhs > sim$h.brk.rhs)
he = sim[he, ]

pe = which(sim$t.brk.lhs < sim$p.brk.lhs | sim$t.brk.rhs > sim$p.brk.rhs)
pe = sim[pe, ]


for (i in 1:nrow(ge)) {
	idx = ge$index[i]
	gi0 = ge$g0[i]
	gi1 = ge$g1[i]
	
	brk.lhs = match(min(ge$t.brk.lhs[i], ge$g.brk.lhs[i]), POS)
	brk.rhs = match(max(ge$t.brk.rhs[i], ge$g.brk.rhs[i]), POS)
	
	g0 = G[, gi0]
	g1 = G[, gi1]
	break
	s0 = g0[brk.lhs : brk.rhs]
	s1 = g1[brk.lhs : brk.rhs]
}



