#
# time & freq dependend genotype pair emissions (LOCAL)
#

library(data.table)


# args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" # args[1] #
tru.file = "../result.truth.local.RData" # args[2] #

sample.n = 10000

load(sprintf("../data.%s.RData", prefix))
load(sprintf("../data.%s.G.RData", prefix)) # genotypes

load(tru.file)



cut.freq = function(f, brk = (0:500)/500, lab = sprintf("%.3f", (1:500)/500)) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

cut.geno = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}

make.matrix = function(len = nrow(G)) {
	matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
}


count.pairs = function(segment, ibd.pair, non.pair, mask = c("00"=1, "01"=2, "02"=3, "11"=4, "12"=5, "22"=6)) {
	ibd = make.matrix()
	non = make.matrix()

	lhs = segment[1] + 5
	rhs = segment[2] - 5

	rng = lhs : rhs

	coord = matrix(c(rng, mask[ibd.pair[rng]]), ncol = 2, byrow = F)
	ibd[ coord ] = ibd[ coord ] + 1

	coord = matrix(c(rng, mask[non.pair[rng]]), ncol = 2, byrow = F)
	non[ coord ] = non[ coord ] + 1

	list(IBD = ibd, NON = non)
}


###


# key = sprintf("%d %d %d %d", truth$h0, truth$h1, truth$lhs.index, truth$rhs.index)
# del = which(duplicated(key))
# tru = truth[-del, ]

tru = truth[sample(nrow(truth)), ]
rm(truth)
gc()


cat("Sampling observed genotypes per site ...\n")

ibd = make.matrix()
non = make.matrix()

j = 0
for (i in 1:nrow(tru)) {
	if (i %% 50 == 0) cat(sprintf(" %d of %d\n", i, sample.n))

	segment = c(tru$lhs.index[i],
							tru$rhs.index[i])

	if (segment[2] - segment[1] < 100) next

	s0 = tru$g0[i]
	s1 = tru$g1[i]
	s0s1 = cut.geno(G[, s0], G[, s1])

	#if (length(which(tru$g0 == s0 & tru$g1 == s1 & tru$lhs.index < segment[2] & tru$rhs.index > segment[1])) > 1) next

	u0 = sample(1:ncol(G), 1)
	u1 = sample(setdiff(1:ncol(G), u0), 1)
	while(length(which(tru$g0 == u0 & tru$g1 == u1 & tru$lhs.index < segment[2] & tru$rhs.index > segment[1])) > 0) {
		u0 = sample(1:ncol(G), 1)
		u1 = sample(setdiff(1:ncol(G), u0), 1)
	}
	u0u1 = cut.geno(G[, u0], G[, u1])

	tmp = count.pairs(segment, s0s1, u0u1)
	ibd = ibd + tmp$IBD
	non = non + tmp$NON

	j = j + 1
	if (j == sample.n) break
}
cat("\n")

count = list(IBD = ibd, NON = non)


pack = paste(sample(c(letters, LETTERS, as.character(0:9)), 16), collapse = "")

save(count, file = sprintf("result.local_ibd_counts.%s.%s.RData", prefix, pack))





###
stop("DONE")
###


library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20.GenErr_1000G"

load(sprintf("data.%s.RData", prefix))



cut.freq = function(f, brk = (0:5000)/5000, lab = sprintf("%.8f", ((1:5000)/5000) - (0.5/5000))) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

make.matrix = function(len = length(POS)) {
	matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
}


res.files = dir(pattern = sprintf("^result\\.local_ibd_counts\\.%s\\..+\\.RData$", prefix), path = "./local_ibd/", full.names = T)

ibd = make.matrix()
non = make.matrix()
for (res.file in res.files) {
	cat(res.file, "\n")

	load(res.file)

	ibd = ibd + count$IBD
	non = non + count$NON
}




cat("Observed genotypes per frequency bin ...\n")

tmp = as.data.frame(ibd)
tmp = split(tmp, cut.freq(AAF))
tmp = lapply(tmp, function(x) {
	x = colSums(x)
	x = x / sum(x)
	as.data.table(as.list(x))
})
tag = names(tmp)
tmp = rbindlist(tmp)
rownames(tmp) = tag
ibd = tmp

tmp = as.data.frame(non)
tmp = split(tmp, cut.freq(AAF))
tmp = lapply(tmp, function(x) {
	x = colSums(x)
	x = x / sum(x)
	as.data.table(as.list(x))
})
tag = names(tmp)
tmp = rbindlist(tmp)
rownames(tmp) = tag
non = tmp

cat("\n")

freq = list(IBD = ibd, NON = non)


### plotting ...

d = NULL

for (t in names(freq)) {
	f = as.numeric(rownames(freq[[t]]))

	for (p in colnames(freq[[t]])) {
		x = freq[[t]][[p]]

		d = rbind(d, data.table(time = t,
														freq = f,
														pair = p,
														prop = x))
	}
}

del = which(is.na(d$prop))
if (length(del) > 0) {
	d = d[-del, ]
}


p = ((1:500)/500) - (0.5/500)
q = 1 - p

expect.non = rbind(data.table(freq = q, pair = "00", prop = p^4),
									 data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
									 data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
									 data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
									 data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
									 data.table(freq = q, pair = "22", prop = q^4))
expect.non$time = "NON"
expect.non$type = "Expected under NON"

expect.ibd = rbind(data.table(freq = q, pair = "00", prop = p^3),
									 data.table(freq = q, pair = "01", prop = 2*(p^2)*q),
									 data.table(freq = q, pair = "02", prop = 0),
									 data.table(freq = q, pair = "11", prop = ((p^2)*q) + (p*(q^2))),
									 data.table(freq = q, pair = "12", prop = 2*p*(q^2)),
									 data.table(freq = q, pair = "22", prop = q^3))
expect.ibd$time = "IBD"
expect.ibd$type = "Expected under IBD"

e = rbind(expect.non, expect.ibd)



# sim$setup = "Simulated, no error"
# err$setup = "Simulated, with error"
# d = rbind(sim, err)


gg = ggplot(d) +
	facet_grid(.~time) +
	geom_point(aes(x = freq, y = prop, colour = pair), size = 0.25, alpha=0.5) +
	geom_line(data = e, aes(x = freq, y = prop, group = pair), size = 0.75, linetype = "22") +
	geom_line(data = e, aes(x = freq, y = prop, colour = pair), size = 0.75, linetype = "22", alpha = 1/3) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Allele frequency") +
	ylab("Observed genotype pair proportion") 

ggsave(filename = "~/Desktop/genotypepairprops.pdf", plot = gg, width = 10, height = 10)


ggplot(d) +
	facet_wrap(~pair) +
	geom_line(aes(x = freq, y = prop, colour = time))






###
stop("DONE")
###


library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20.GenErr_1000G" # args[1] #

sample.n = 10000

load(sprintf("./data.%s.RData", prefix))
load(sprintf("./data.%s.G.RData", prefix)) # genotypes




make.matrix = function(len = nrow(G)) {
	matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
}


cut.freq = function(f, brk = (0:5000)/5000, lab = sprintf("%.8f", ((1:5000)/5000) - (0.5/5000))) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}


cut.geno = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}


count.pairs = function(n, pair, mask = c("00"=1, "01"=2, "02"=3, "11"=4, "12"=5, "22"=6)) {
	count = make.matrix()

	coord = matrix(c(1:n, mask[pair]), ncol = 2, byrow = F)
	count[ coord ] = count[ coord ] + 1
	
	count
}



cat("Sampling observed genotypes per site ...\n")

count = make.matrix()

j = 0
while(T) {
	if (j %% 100 == 0) cat(sprintf(" %d of %d\n", j, sample.n))
	
	i0 = sample(1:ncol(G), 1)
	i1 = sample(setdiff(1:ncol(G), i0), 1)

	i0i1 = cut.geno(G[, i0], G[, i1])
	
	count = count + count.pairs(length(POS), i0i1)
	
	j = j + 1
	if (j == sample.n) break
}
cat("\n")



tmp = data.table(g0 = GF0, g1 = GF1, g2 = GF2, e0 = AAF^2, e1 = 2 *  AAF * (1-AAF), e1 = (1-AAF)^2)

count = data.table(g00 = GF0^2, g01 = 2 * GF0 * GF1, g02 = 2 * GF0 * GF2, g11 = 4 * GF0 * GF2, g12 = 2 * GF2 * GF1, g22 = GF2^2)

p = 1-AAF
q = AAF

tmp = data.table(g00 = p^4, g01 = 4 * p^3 * q, g02 = 2 * GF0 * GF2, g11 = 4 * GF0 * GF2, g12 = 2 * GF2 * GF1, g22 = GF2^2)


tmp = as.data.frame(count)
cut = cut.freq(AAF)
tmp = split(tmp, cut)
tmp = lapply(tmp, function(x) {
  x = colSums(x)
	x = x / sum(x)
	as.data.table(as.list(x))
})
tmp = rbindlist(tmp)

d = NULL
f = as.numeric(levels(cut))
for (p in names(tmp)) {
	x = tmp[[p]]
	
	d = rbind(d, data.table(freq = f,
													pair = p,
													prop = x))
}

ggplot(d) +	geom_point(aes(x = freq, y = prop, colour = pair), size = 0.25, alpha=0.5)



