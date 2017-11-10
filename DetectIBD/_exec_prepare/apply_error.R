#
# apply error profile
#

library(ggplot2)
library(ggthemes)
library(data.table)

#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]
err.file = "_exec_prepare/result.profile.1000g.platinum.RData" #args[2]
out.name = "GenErr_1000G" #args[3]


load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.H.RData", prefix)) # original haplotypes
load(sprintf("data.%s.G.RData", prefix)) # original genotypes
load(err.file)



dir.create(out.name, showWarnings = F)
setwd(out.name)
cat("Data in: ", getwd(), "\n")




###
approx.error = function(af, tab = error.table.empirical) {
	a00 = approx(tab$freq, tab$p00, xout = af, rule = 2)
	a01 = approx(tab$freq, tab$p01, xout = af, rule = 2)
	a02 = approx(tab$freq, tab$p02, xout = af, rule = 2)

	a10 = approx(tab$freq, tab$p10, xout = af, rule = 2)
	a11 = approx(tab$freq, tab$p11, xout = af, rule = 2)
	a12 = approx(tab$freq, tab$p12, xout = af, rule = 2)

	a20 = approx(tab$freq, tab$p20, xout = af, rule = 2)
	a21 = approx(tab$freq, tab$p21, xout = af, rule = 2)
	a22 = approx(tab$freq, tab$p22, xout = af, rule = 2)

	n0 = a00$y + a01$y + a02$y
	n1 = a10$y + a11$y + a12$y
	n2 = a20$y + a21$y + a22$y

	data.table(freq = af,
						 p00 = a00$y / n0, p01 = a01$y / n0, p02 = a02$y / n0,
						 p10 = a10$y / n1, p11 = a11$y / n1, p12 = a12$y / n1,
						 p20 = a20$y / n2, p21 = a21$y / n2, p22 = a22$y / n2)
}




### apply error to haplotypes
cat("Apply error profile to haplotypes\n")

idx.m = (1:ncol(H))[(1:ncol(H)) %% 2 != 0]
idx.p = (1:ncol(H))[(1:ncol(H)) %% 2 == 0]

tpl.m = c(0, 0, 1, 1)
tpl.p = c(0, 1, 0, 1)

for (i in 1:nrow(H)) {
	if (i %% 5000 == 0) cat(sprintf(" %d of %d\n", i, nrow(H), i / nrow(H) * 100))

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

	err = approx.error(AAF[i])
	if (any(is.na(err))) {
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


save(H, file = sprintf("data.%s.%s.H.RData", prefix, out.name))



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
	if (i %% 5000 == 0) cat(sprintf(" %d of %d\n", i, nrow(H), i / nrow(H) * 100))

	g = G[i, ]

	g0 = which(g == 0)
	g1 = which(g == 1)
	g2 = which(g == 2)

	h = H[i, ]
	q = h[idx.m] + h[idx.p]

	if (length(g0) > 0) {
		x = tabulate(q[g0] + 1, nbins = 3)
		x00[i] = x[1]
		x01[i] = x[2]
		x02[i] = x[3]
	}
	if (length(g1) > 0) {
		x = tabulate(q[g1] + 1, nbins = 3)
		x10[i] = x[1]
		x11[i] = x[2]
		x12[i] = x[3]
	}
	if (length(g2) > 0) {
		x = tabulate(q[g2] + 1, nbins = 3)
		x20[i] = x[1]
		x21[i] = x[2]
		x22[i] = x[3]
	}
}

applied.error.table = data.table(freq = AAF,
																 e00 = x00, e01 = x01, e02 = x02,
																 e10 = x10, e11 = x11, e12 = x12,
																 e20 = x20, e21 = x21, e22 = x22)
cat("OK\n")


ae = split(applied.error.table, applied.error.table$freq)
ae = lapply(ae, function(ae) {
	data.table(freq = ae$freq[1],
						 e00 = sum(ae$e00), e01 = sum(ae$e01), e02 = sum(ae$e02),
						 e10 = sum(ae$e10), e11 = sum(ae$e11), e12 = sum(ae$e12),
						 e20 = sum(ae$e20), e21 = sum(ae$e21), e22 = sum(ae$e22))
})
ae = rbindlist(ae)

n0 = ae$e00 + ae$e01 + ae$e02
n1 = ae$e10 + ae$e11 + ae$e12
n2 = ae$e20 + ae$e21 + ae$e22

ae = rbind(
	data.table(freq = ae$freq, true.gt = "0", call.gt = "0", prop = ae$e00 / n0),
	data.table(freq = ae$freq, true.gt = "0", call.gt = "1", prop = ae$e01 / n0),
	data.table(freq = ae$freq, true.gt = "0", call.gt = "2", prop = ae$e02 / n0),

	data.table(freq = ae$freq, true.gt = "1", call.gt = "0", prop = ae$e10 / n1),
	data.table(freq = ae$freq, true.gt = "1", call.gt = "1", prop = ae$e11 / n1),
	data.table(freq = ae$freq, true.gt = "1", call.gt = "2", prop = ae$e12 / n1),

	data.table(freq = ae$freq, true.gt = "2", call.gt = "0", prop = ae$e20 / n2),
	data.table(freq = ae$freq, true.gt = "2", call.gt = "1", prop = ae$e21 / n2),
	data.table(freq = ae$freq, true.gt = "2", call.gt = "2", prop = ae$e22 / n2)
)


ee = approx.error(sort(unique(AAF)))
ee = approx.error((0:5000) / 5000)

n0 = ee$p00 + ee$p01 + ee$p02
n1 = ee$p10 + ee$p11 + ee$p12
n2 = ee$p20 + ee$p21 + ee$p22

ee = rbind(
	data.table(freq = ee$freq, true.gt = "0", call.gt = "0", prop = ee$p00 / n0),
	data.table(freq = ee$freq, true.gt = "0", call.gt = "1", prop = ee$p01 / n0),
	data.table(freq = ee$freq, true.gt = "0", call.gt = "2", prop = ee$p02 / n0),

	data.table(freq = ee$freq, true.gt = "1", call.gt = "0", prop = ee$p10 / n1),
	data.table(freq = ee$freq, true.gt = "1", call.gt = "1", prop = ee$p11 / n1),
	data.table(freq = ee$freq, true.gt = "1", call.gt = "2", prop = ee$p12 / n1),

	data.table(freq = ee$freq, true.gt = "2", call.gt = "0", prop = ee$p20 / n2),
	data.table(freq = ee$freq, true.gt = "2", call.gt = "1", prop = ee$p21 / n2),
	data.table(freq = ee$freq, true.gt = "2", call.gt = "2", prop = ee$p22 / n2)
)

ee$true.gt = sprintf("True genotype = %s", ee$true.gt)
ee$call.gt = sprintf("Observed genotype = %s", ee$call.gt)

ae$true.gt = sprintf("True genotype = %s", ae$true.gt)
ae$call.gt = sprintf("Observed genotype = %s", ae$call.gt)

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


gg = ggplot(data = ee) +
	facet_wrap(~true.gt, nrow = 1) +
	geom_point(aes(x = freq, y = prop, colour=call.gt), size = 1, shape=16, alpha=1/3) +
	scale_x_continuous(expand = c(0.025, 0)) +
	scale_y_continuous(expand = c(0.025, 0)) +
	scale_colour_colorblind() +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.margin.x = unit(x = 0.25, units = "inches"),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Allele frequency") +
	ylab("Relative error proportion")
ggsave(gg, filename = sprintf("_plot.errorprofile.%s.pdf", out.name), width = 12, height = 10)
ggsave(gg, filename = sprintf("_plot.errorprofile.%s.png", out.name), width = 12, height = 10)



# make genotypes
# keep track of introduced error
cat("Make genotypes from error-treated haplotypes\n")
applied.error.matrix = matrix(0, nrow = 3, ncol = 3, dimnames = list(as.character(0:2), as.character(0:2)))
j = 1
for (i in 1:ncol(G)) {
	if (i %% 100 == 0) cat(".")

	g0 = which(G[, i] == 0)
	g1 = which(G[, i] == 1)
	g2 = which(G[, i] == 2)

	G[, i] = H[, j] + H[, j+1]
	j = j + 2

	applied.error.matrix["0", "0"] = applied.error.matrix["0", "0"] + sum(G[g0, i] == 0)
	applied.error.matrix["0", "1"] = applied.error.matrix["0", "1"] + sum(G[g0, i] == 1)
	applied.error.matrix["0", "2"] = applied.error.matrix["0", "2"] + sum(G[g0, i] == 2)

	applied.error.matrix["1", "0"] = applied.error.matrix["1", "0"] + sum(G[g1, i] == 0)
	applied.error.matrix["1", "1"] = applied.error.matrix["1", "1"] + sum(G[g1, i] == 1)
	applied.error.matrix["1", "2"] = applied.error.matrix["1", "2"] + sum(G[g1, i] == 2)

	applied.error.matrix["2", "0"] = applied.error.matrix["2", "0"] + sum(G[g2, i] == 0)
	applied.error.matrix["2", "1"] = applied.error.matrix["2", "1"] + sum(G[g2, i] == 1)
	applied.error.matrix["2", "2"] = applied.error.matrix["2", "2"] + sum(G[g2, i] == 2)
}
cat("\n")


save(G, file = sprintf("data.%s.%s.G.RData", prefix, out.name))


em = error.matrix / rowSums(error.matrix) * 100
em = rbind(data.table(true.gt = "0", call.gt = "0", txt = sprintf("%.4f%%", em["0", "0"])),
					 data.table(true.gt = "0", call.gt = "1", txt = sprintf("%.4f%%", em["0", "1"])),
					 data.table(true.gt = "0", call.gt = "2", txt = sprintf("%.4f%%", em["0", "2"])),

					 data.table(true.gt = "1", call.gt = "0", txt = sprintf("%.4f%%", em["1", "0"])),
					 data.table(true.gt = "1", call.gt = "1", txt = sprintf("%.4f%%", em["1", "1"])),
					 data.table(true.gt = "1", call.gt = "2", txt = sprintf("%.4f%%", em["1", "2"])),

					 data.table(true.gt = "2", call.gt = "0", txt = sprintf("%.4f%%", em["2", "0"])),
					 data.table(true.gt = "2", call.gt = "1", txt = sprintf("%.4f%%", em["2", "1"])),
					 data.table(true.gt = "2", call.gt = "2", txt = sprintf("%.4f%%", em["2", "2"])))

am = applied.error.matrix / rowSums(applied.error.matrix) * 100
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
MAC = AAC = rowSums(H)
MAF = AAF = AAC / ncol(H)
i = which(MAC > (ncol(H) / 2))
if (length(i) > 0) {
	MAC[i] = ncol(H) - MAC[i]
	MAF[i] = 1 - MAF[i]
}

GC0 = apply(G, 1, function(x) length(which(x == 0)))
GC1 = apply(G, 1, function(x) length(which(x == 1)))
GC2 = apply(G, 1, function(x) length(which(x == 2)))

GF0 = GC0 / ncol(G)
GF1 = GC1 / ncol(G)
GF2 = GC2 / ncol(G)



# make fk index
cat("Compiling rare variant index\n")
fklist = 2:25 # c(2:25, seq(30, 45, by = 5), seq(50, 90, by = 10), seq(100, 500, by=50))
FKI = list()
cum = 0
i = 0
for (fk in fklist) {
	i = i + 1
	idx = which(AAC == fk)
	num = length(idx)
	cmb = choose(fk, 2) * num
	cum = cmb + cum
	cat(sprintf("# f_%d : %d (%d pairs | %d cumulative pairs)\n", fk, num, cmb, cum))

	tmp = lapply(idx, function(x) {
		g = which(G[x, ] == 1)
		if (length(g) < 2) return(NULL) # others may be autozygous, hence higher count
		as.data.table(cbind(x, POS[x], AAC[x], t(combn(g, 2))))
	})
	tmp = rbindlist(tmp)
	names(tmp) = c("index", "pos", "fk", "g0", "g1")

	a = tmp$g0 * 2 - 1
	b = tmp$g0 * 2
	coord.a = as.matrix(data.frame(x = tmp$index, y = a))
	coord.b = as.matrix(data.frame(x = tmp$index, y = b))
	coord.a = H[ coord.a ]
	coord.b = H[ coord.b ]
	if (any(coord.a == coord.b)) stop("???")
	tmp$h0 = (a * coord.a) + (b * coord.b)

	a = tmp$g1 * 2 - 1
	b = tmp$g1 * 2
	coord.a = as.matrix(data.frame(x = tmp$index, y = a))
	coord.b = as.matrix(data.frame(x = tmp$index, y = b))
	coord.a = H[ coord.a ]
	coord.b = H[ coord.b ]
	if (any(coord.a == coord.b)) stop("???")
	tmp$h1 = (a * coord.a) + (b * coord.b)

	FKI = rbind(FKI, tmp)
}



# save data
cat("Saving data\n")
save(POS, CPOS, BRK, CBRK,
		 AAC, AAF, MAC, MAF,
		 GC0, GC1, GC2, GF0, GF1, GF2,
		 FKI,
		 applied.error.table, applied.error.matrix,
		 file = sprintf("data.%s.%s.RData", prefix, out.name))

cat("OK\n")






