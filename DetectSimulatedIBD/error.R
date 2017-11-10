#
# apply error profile
#

library(ggplot2)
library(ggthemes)
library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]
err.file = args[2]
err.data = args[3]
err.func = args[4] ### 'approx' or 'model'
out.name = args[5]


load(sprintf("%s.RData", prefix))
load(err.file)


error.data = NULL
error.func = function() {}
error.matx = get(sprintf("err.mat.%s", err.data))

if (err.func == "approx") {
	error.data = get(sprintf("err.frq.%s", err.data))
	error.func = function(data, af, gf0, gf1, gf2) {
		approx.error(data, af, plotting = F)
	}
} else {
	if (err.func == "model") {
		error.data = get(sprintf("err.mat.%s", err.data))
		error.func = function(data, af, gf0, gf1, gf2) {
			model.error(data, gf0, gf1, gf2, plotting = F)
		}
	} else {
		stop("Unknown data/function")
	}
}


H = as.matrix(load.bigmatrix(sprintf("%s.H", prefix)))
G = as.matrix(load.bigmatrix(sprintf("%s.G", prefix)))

AC = rowSums(H)
AF = AC / ncol(H)

gf0 = apply(G, 1, function(x) length(which(x == 0)) / 2500)
gf1 = apply(G, 1, function(x) length(which(x == 1)) / 2500)
gf2 = apply(G, 1, function(x) length(which(x == 2)) / 2500)


dir.create(out.name, showWarnings = F)
setwd(out.name)
cat("Data in: ", getwd(), "\n")



### apply error to haplotypes
cat("Apply error profile to haplotypes\n")

idx.m = (1:ncol(H))[(1:ncol(H)) %% 2 != 0]
idx.p = (1:ncol(H))[(1:ncol(H)) %% 2 == 0]

tpl.m = c(0, 0, 1, 1)
tpl.p = c(0, 1, 0, 1)

for (i in 1:nrow(H)) {
	if (i %% 10000 == 0) cat(".")
	
	h = H[i, ]
	hm = h[idx.m]
	hp = h[idx.p]
	
	i00 = which(hm == 0 & hp == 0)
	i01 = which(hm == 0 & hp == 1)
	i10 = which(hm == 1 & hp == 0)
	i11 = which(hm == 1 & hp == 1)
	
	n00 = length(i00)
	n01 = length(i01)
	n10 = length(i10)
	n11 = length(i11)
	sum = n00 + n01 + n10 + n11
	
	err = error.func(error.data, AF[i], n00 / sum, (n01 + n10) / sum, n11 / sum)
	if (nrow(err) == 0) {
		next
	}
	
	if (n00 != 0) {
		x = sample(c(1, 2, 3, 4), n00, replace = T, prob = c(err$p00, err$p01 / 2, err$p01 / 2, err$p02))
		hm[i00] = tpl.m[x]
		hp[i00] = tpl.p[x]
	}
	if (n01 != 0) {
		x = sample(c(1, 2,    4), n01, replace = T, prob = c(err$p10, err$p11    ,              err$p12))
		hm[i01] = tpl.m[x]
		hp[i01] = tpl.p[x]
	}
	if (n10 != 0) {
		x = sample(c(1,    3, 4), n10, replace = T, prob = c(err$p10             , err$p11    , err$p12))
		hm[i10] = tpl.m[x]
		hp[i10] = tpl.p[x]
	}
	if (n11 != 0) {
		x = sample(c(1, 2, 3, 4), n11, replace = T, prob = c(err$p20, err$p21 / 2, err$p21 / 2, err$p22))
		hm[i11] = tpl.m[x]
		hp[i11] = tpl.p[x]
	}
	
	h[idx.m] = hm
	h[idx.p] = hp
	H[i, ] = h
}
cat("OK\n")


# measure applied error
cat("Measure error\n")

idx.m = (1:ncol(H))[(1:ncol(H)) %% 2 != 0]
idx.p = (1:ncol(H))[(1:ncol(H)) %% 2 == 0]

x00 = rep(0, nrow(H))
x01 = rep(0, nrow(H))
x02 = rep(0, nrow(H))

x10 = rep(0, nrow(H))
x11 = rep(0, nrow(H))
x12 = rep(0, nrow(H))

x20 = rep(0, nrow(H))
x21 = rep(0, nrow(H))
x22 = rep(0, nrow(H))

for (i in 1:nrow(H)) {
	if (i %% 10000 == 0) cat(".")
	
	g = G[i, ]
	
	g0 = which(g == 0)
	g1 = which(g == 1)
	g2 = which(g == 2)
	
	h = H[i, ]
	h = h[idx.m] + h[idx.p]
	
	if (length(g0) > 0) {
		x = tabulate(h[g0] + 1, nbins = 3)
		x00[i] = x[1]
		x01[i] = x[2]
		x02[i] = x[3]
	}
	if (length(g1) > 0) {
		x = tabulate(h[g1] + 1, nbins = 3)
		x10[i] = x[1]
		x11[i] = x[2]
		x12[i] = x[3]
	}
	if (length(g2) > 0) {
		x = tabulate(h[g2] + 1, nbins = 3)
		x20[i] = x[1]
		x21[i] = x[2]
		x22[i] = x[3]
	}
}

applied.err.frq = data.table(af = AF,
														 from.0.to.0 = x00, from.0.to.1 = x01, from.0.to.2 = x02,
														 from.1.to.0 = x10, from.1.to.1 = x11, from.1.to.2 = x12,
														 from.2.to.0 = x20, from.2.to.1 = x21, from.2.to.2 = x22)
cat("OK\n")


ae = split(applied.err.frq, applied.err.frq$af)
ae = lapply(ae, function(ae) {
	data.table(af = ae$af[1],
						 from.0.to.0 = sum(ae$from.0.to.0), from.0.to.1 = sum(ae$from.0.to.1), from.0.to.2 = sum(ae$from.0.to.2),
						 from.1.to.0 = sum(ae$from.1.to.0), from.1.to.1 = sum(ae$from.1.to.1), from.1.to.2 = sum(ae$from.1.to.2),
						 from.2.to.0 = sum(ae$from.2.to.0), from.2.to.1 = sum(ae$from.2.to.1), from.2.to.2 = sum(ae$from.2.to.2))
})
ae = rbindlist(ae)

n0 = ae$from.0.to.0 + ae$from.0.to.1 + ae$from.0.to.2
n1 = ae$from.1.to.0 + ae$from.1.to.1 + ae$from.1.to.2
n2 = ae$from.2.to.0 + ae$from.2.to.1 + ae$from.2.to.2
ae = rbind(
	data.table(freq = ae$af, true.gt = "0", call.gt = "0", prop = ae$from.0.to.0 / n0),
	data.table(freq = ae$af, true.gt = "0", call.gt = "1", prop = ae$from.0.to.1 / n0),
	data.table(freq = ae$af, true.gt = "0", call.gt = "2", prop = ae$from.0.to.2 / n0),
	
	data.table(freq = ae$af, true.gt = "1", call.gt = "0", prop = ae$from.1.to.0 / n1),
	data.table(freq = ae$af, true.gt = "1", call.gt = "1", prop = ae$from.1.to.1 / n1),
	data.table(freq = ae$af, true.gt = "1", call.gt = "2", prop = ae$from.1.to.2 / n1),
	
	data.table(freq = ae$af, true.gt = "2", call.gt = "0", prop = ae$from.2.to.0 / n2),
	data.table(freq = ae$af, true.gt = "2", call.gt = "1", prop = ae$from.2.to.1 / n2),
	data.table(freq = ae$af, true.gt = "2", call.gt = "2", prop = ae$from.2.to.2 / n2)
)

ee = error.func(error.data, AF, gf0, gf1, gf2)
ee = split(ee, ee$freq)
ee = lapply(ee, function(ee) {
	data.table(af = ee$freq[1],
						 p00 = sum(ee$p00), p01 = sum(ee$p01), p02 = sum(ee$p02),
						 p10 = sum(ee$p10), p11 = sum(ee$p11), p12 = sum(ee$p12),
						 p20 = sum(ee$p20), p21 = sum(ee$p21), p22 = sum(ee$p22))
})
ee = rbindlist(ee)

n0 = ee$p00 + ee$p01 + ee$p02
n1 = ee$p10 + ee$p11 + ee$p12
n2 = ee$p20 + ee$p21 + ee$p22
ee = rbind(
	data.table(freq = ee$af, true.gt = "0", call.gt = "0", prop = ee$p00 / n0),
	data.table(freq = ee$af, true.gt = "0", call.gt = "1", prop = ee$p01 / n0),
	data.table(freq = ee$af, true.gt = "0", call.gt = "2", prop = ee$p02 / n0),
	
	data.table(freq = ee$af, true.gt = "1", call.gt = "0", prop = ee$p10 / n1),
	data.table(freq = ee$af, true.gt = "1", call.gt = "1", prop = ee$p11 / n1),
	data.table(freq = ee$af, true.gt = "1", call.gt = "2", prop = ee$p12 / n1),
	
	data.table(freq = ee$af, true.gt = "2", call.gt = "0", prop = ee$p20 / n2),
	data.table(freq = ee$af, true.gt = "2", call.gt = "1", prop = ee$p21 / n2),
	data.table(freq = ee$af, true.gt = "2", call.gt = "2", prop = ee$p22 / n2)
)

ee$true.gt = sprintf("True genotype is %s", ee$true.gt)
ee$call.gt = sprintf("Observed genotype is %s", ee$call.gt)

ae$true.gt = sprintf("True genotype is %s", ae$true.gt)
ae$call.gt = sprintf("Observed genotype is %s", ae$call.gt)

gg = ggplot(data = ee) + 
	facet_grid(call.gt~true.gt) + 
	geom_point(data = ae, aes(x = freq, y = prop), size = 1.25, shape=16, alpha=1/3, colour="orangered") +
	geom_point(data = ee, aes(x = freq, y = prop), size = 0.75, shape=16, alpha=1/3, colour="blue") + 
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Allele frequency") +
	ylab("Relative error proportion")



# make genotypes
# keep track of introduced error
cat("Make genotypes from errorised haplotypes\n")
applied.err.mat = matrix(0, nrow = 3, ncol = 3, dimnames = list(as.character(0:2), as.character(0:2)))
j = 1
for (i in 1:ncol(G)) {
	if (i %% 100 == 0) cat(".")
	
	g0 = which(G[, i] == 0)
	g1 = which(G[, i] == 1)
	g2 = which(G[, i] == 2)
	
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
	
	applied.err.mat["0", "0"] = applied.err.mat["0", "0"] + sum(G[g0, i] == 0)
	applied.err.mat["0", "1"] = applied.err.mat["0", "1"] + sum(G[g0, i] == 1)
	applied.err.mat["0", "2"] = applied.err.mat["0", "2"] + sum(G[g0, i] == 2)
	
	applied.err.mat["1", "0"] = applied.err.mat["1", "0"] + sum(G[g1, i] == 0)
	applied.err.mat["1", "1"] = applied.err.mat["1", "1"] + sum(G[g1, i] == 1)
	applied.err.mat["1", "2"] = applied.err.mat["1", "2"] + sum(G[g1, i] == 2)
	
	applied.err.mat["2", "0"] = applied.err.mat["2", "0"] + sum(G[g2, i] == 0)
	applied.err.mat["2", "1"] = applied.err.mat["2", "1"] + sum(G[g2, i] == 1)
	applied.err.mat["2", "2"] = applied.err.mat["2", "2"] + sum(G[g2, i] == 2)
}
cat("\n")




em = error.matx / rowSums(error.matx) * 100
em = rbind(data.table(true.gt = "0", call.gt = "0", txt = sprintf("%.4f%%", em["0", "0"])),
					 data.table(true.gt = "0", call.gt = "1", txt = sprintf("%.4f%%", em["0", "1"])),
					 data.table(true.gt = "0", call.gt = "2", txt = sprintf("%.4f%%", em["0", "2"])),
					 
					 data.table(true.gt = "1", call.gt = "0", txt = sprintf("%.4f%%", em["1", "0"])),
					 data.table(true.gt = "1", call.gt = "1", txt = sprintf("%.4f%%", em["1", "1"])),
					 data.table(true.gt = "1", call.gt = "2", txt = sprintf("%.4f%%", em["1", "2"])),
					 
					 data.table(true.gt = "2", call.gt = "0", txt = sprintf("%.4f%%", em["2", "0"])),
					 data.table(true.gt = "2", call.gt = "1", txt = sprintf("%.4f%%", em["2", "1"])),
					 data.table(true.gt = "2", call.gt = "2", txt = sprintf("%.4f%%", em["2", "2"])))

am = applied.err.mat / rowSums(applied.err.mat) * 100
am = rbind(data.table(true.gt = "0", call.gt = "0", txt = sprintf("%.4f%%", am["0", "0"])),
					 data.table(true.gt = "0", call.gt = "1", txt = sprintf("%.4f%%", am["0", "1"])),
					 data.table(true.gt = "0", call.gt = "2", txt = sprintf("%.4f%%", am["0", "2"])),
					 
					 data.table(true.gt = "1", call.gt = "0", txt = sprintf("%.4f%%", am["1", "0"])),
					 data.table(true.gt = "1", call.gt = "1", txt = sprintf("%.4f%%", am["1", "1"])),
					 data.table(true.gt = "1", call.gt = "2", txt = sprintf("%.4f%%", am["1", "2"])),
					 
					 data.table(true.gt = "2", call.gt = "0", txt = sprintf("%.4f%%", am["2", "0"])),
					 data.table(true.gt = "2", call.gt = "1", txt = sprintf("%.4f%%", am["2", "1"])),
					 data.table(true.gt = "2", call.gt = "2", txt = sprintf("%.4f%%", am["2", "2"])))

em$true.gt = sprintf("True genotype is %s", em$true.gt)
em$call.gt = sprintf("Observed genotype is %s", em$call.gt)

am$true.gt = sprintf("True genotype is %s", am$true.gt)
am$call.gt = sprintf("Observed genotype is %s", am$call.gt)

gg = gg + geom_label(data = em, aes(x = 0.5, y = 0.54, label = txt), size = 3.5, colour = "blue")
gg = gg + geom_label(data = am, aes(x = 0.5, y = 0.46, label = txt), size = 3.5, colour = "orangered")

ggsave(filename = sprintf("_plot.errorprofile.%s.%s.pdf", prefix, out.name), plot = gg, width = 10, height = 10)


# get frequencies
cat("Getting frequencies\n")
MAC = rowSums(H)
MAF = MAC / ncol(H)
FRQ = MAF
i = which(MAC > (ncol(H) / 2))
if (length(i) > 0) {
	MAC[i] = ncol(H) - MAC[i]
	MAF[i] = 1 - MAF[i]
}


# recover previouslty defined fk max
fk.max = max(FKI$n.sharer)

# make fk index
cat("Compiling rare variant index\n")
FKI = which(MAC > 1 & MAC <= fk.max)
cat("# rare variants: ", length(FKI))
FKI = lapply(FKI, function(i) {
	g = which(G[i, ] == 1)
	h = if (FRQ[i] <= 0.5) which(H[i, ] == 1) else which(H[i, ] == 0)
	h = intersect(h, convert.gen2hap.sharers(g))
	n = length(g)
	if (n < 2) return(NULL)
	data.table(index = i, 
						 pos = POS[i], 
						 frq = FRQ[i],
						 maf = MAF[i],
						 mac = MAC[i],
						 n.sharer = n, 
						 g.sharer = paste(as.character(g), collapse = "|"),
						 h.sharer = paste(as.character(h), collapse = "|"))
})
FKI = rbindlist(FKI)



# storing on disk
cat("Storing data on disk\n")

if (file.exists(sprintf("_bigmatrix.%s.%s.H.bin",  prefix, out.name)) ||
		file.exists(sprintf("_bigmatrix.%s.%s.H.desc", prefix, out.name))) {
	unlink(c(sprintf("_bigmatrix.%s.%s.H.bin",  prefix, out.name), 
					 sprintf("_bigmatrix.%s.%s.H.desc", prefix, out.name)))
}

tmp = big.matrix(nrow(H), ncol(H), type = "char", 
								 backingpath = getwd(),
								 backingfile =    sprintf("_bigmatrix.%s.%s.H.bin",  prefix, out.name), 
								 descriptorfile = sprintf("_bigmatrix.%s.%s.H.desc", prefix, out.name))

for (i in 1:ncol(H)) {
	tmp[, i] = H[, i]
	if (i %% 100 == 0) cat(".")
}

if (file.exists(sprintf("_bigmatrix.%s.%s.G.bin",  prefix, out.name)) ||
		file.exists(sprintf("_bigmatrix.%s.%s.G.desc", prefix, out.name))) {
	unlink(c(sprintf("_bigmatrix.%s.%s.G.bin",  prefix, out.name), 
					 sprintf("_bigmatrix.%s.%s.G.desc", prefix, out.name)))
}

tmp = big.matrix(nrow(G), ncol(G), type = "char", 
								 backingpath = getwd(),
								 backingfile =    sprintf("_bigmatrix.%s.%s.G.bin",  prefix, out.name), 
								 descriptorfile = sprintf("_bigmatrix.%s.%s.G.desc", prefix, out.name))

for (i in 1:ncol(G)) {
	tmp[, i] = G[, i]
	if (i %% 100 == 0) cat(".")
}

H = sprintf("%s.%s.H", prefix, out.name)
G = sprintf("%s.%s.G", prefix, out.name)


P = sprintf("%s.%s.P", prefix, out.name)


# save data
cat("Saving data\n")
save(H, G, P,
		 POS, MAC, MAF, FKI, 
		 parse.sharers,
		 convert.gen2hap.sharers,
		 load.bigmatrix,
		 applied.err.mat, applied.err.frq, err.data, err.func, ### <--
		 file = sprintf("%s.%s.RData", prefix, out.name))




