#
# plot detected vs truth
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


prefix = "beagle_ibd.history.generror_affymetrixaxiom.P"

beagle.file = sprintf("%s.ibd", prefix)
sim.file = "history.generror_affymetrixaxiom.RData"

truth.file = "../truth/ibd_mrca.truth.RData"


load(truth.file)
load(sim.file)


pos = round(POS)

i = anyDuplicated(pos)
while (i != 0) {
	if (i == 1) {
		pos[i] = pos[i] - 1
	} else {
		pos[i] = pos[i] + 1
	}
	i = anyDuplicated(pos)
}

truth$position = pos[match(truth$position, POS)]
truth$lhs = pos[match(truth$lhs, POS)]
truth$rhs = pos[match(truth$rhs, POS)]


ibd = read.table(beagle.file, header = F, stringsAsFactors = F)
names(ibd) = c("indv1", "x1", "indv2", "x2", "chr", "beg", "end", "lod")

ibd$g0 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd$indv1))
ibd$g1 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd$indv2))

ibd$h0 = (ibd$g0 * 2) + (ibd$x1 - 2)
ibd$h1 = (ibd$g1 * 2) + (ibd$x2 - 2)


ibd.key = sprintf("%d %d", ibd$h0, ibd$h1)
tru.key = sprintf("%d %d", truth$h0, truth$h1)
key = intersect(ibd.key, tru.key)

ibd = ibd[match(key, ibd.key), ]
tru = truth[match(key, tru.key), ]


a = which(tru$position >= ibd$beg)
b = which(tru$position <= ibd$end)
x = intersect(a, b)

match = data.table(position = tru$position[x],
									 fk = tru$fk[x],
									 wall.lhs = tru$wall.lhs[x],
									 wall.rhs = tru$wall.rhs[x],
									 t.lhs = tru$lhs[x],
									 t.rhs = tru$rhs[x],
									 b.lhs = ibd$beg[x],
									 b.rhs = ibd$end[x])


del = which(match$wall.lhs == T); if (length(del) > 0) match = match[-del, ]
del = which(match$wall.rhs == T); if (length(del) > 0) match = match[-del, ]

del = which(match$b.lhs < 100); if (length(del) > 0) match = match[-del, ]
del = which(match$b.rhs > max(pos) - 100); if (length(del) > 0) match = match[-del, ]

match$t.brk.lhs = pos[ match(match$t.lhs, pos) - 1 ]
match$t.brk.rhs = pos[ match(match$t.rhs, pos) + 1 ]
match$b.brk.lhs = pos[ match(match$b.lhs, pos) - 1 ]
match$b.brk.rhs = pos[ match(match$b.rhs, pos) + 1 ]


match$t.brk.dist.lhs = match$t.brk.lhs - match$position
match$t.brk.dist.rhs = match$t.brk.rhs - match$position

match$b.brk.dist.lhs = match$b.brk.lhs - match$position
match$b.brk.dist.rhs = match$b.brk.rhs - match$position



d = rbind(data.table(side="LHS", fk = match$fk, det = match$b.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, det = match$b.brk.dist.rhs, tru = match$t.brk.dist.rhs))



### mapped

d$map = (d$det / d$tru)
z = which(d$side == "LHS")
d$map[z] = d$map[z] * -1


lhs = density(d$map, bw = 0.005, n = 1024*10, from=-5, to=0)
lhs = data.table(x = lhs$x, y = lhs$y / sum(lhs$y))
rhs = density(d$map, bw = 0.005, n = 1024*10, from=0, to=5)
rhs = data.table(x = rhs$x, y = rhs$y / sum(rhs$y))
p = rbind(lhs, rhs)


gg = ggplot(data=p) + 
	geom_vline(xintercept=c(-1, 0, 1), colour="grey") +
	geom_line(aes(x=x, y=y), alpha=0.75, size=1) + 
	scale_x_continuous(breaks = -5:5) +
	theme_few() +
	theme(aspect.ratio=1/3,
				legend.title=element_blank(),
				legend.position = "top") +
	xlab("Relative physical distance") +
	ylab("Density")

ggsave(filename = sprintf("_plot.beagle_ibd-truth.mapped.%s.pdf", prefix), plot = gg, width = 15, height = 5)


### scatter

brk = seq(log10(1), log10(100e6), length.out = 251)

p = lapply(split(d, d$side), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$side  = x$side[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


p$Var1 = p$Var1 + 1
x = (p$side == "LHS")
p$Var1[x] = p$Var1[x] * -1


tag = c("100b", "10Kb", "1Mb", "100Mb")
tik = as.numeric(cut(log10(c(100, 10000, 1e6, 100e6)), breaks = brk, include.lowest = T)) + 1


gg = ggplot(data = p) + 
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(-1, 1), intercept = c(-2, -2), colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik), tik), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1/2) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance (IBD)")

ggsave(filename = sprintf("_plot.beagle_ibd-truth.scatter.%s.pdf", prefix), plot = gg, width = 15, height = 5)




