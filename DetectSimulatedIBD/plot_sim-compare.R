#
# plot: compare sim/err data
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)

args = commandArgs(T)

sim.history = args[1] 
err.history = args[2] 
compare.file = args[3]


curdir = getwd()
setwd(dirname(sim.history))

sim.prefix = sub("^(.+)\\.RData$", "\\1", basename(sim.history))
load(sprintf("%s.RData", sim.prefix))

SG = load.bigmatrix(sprintf("%s.G", sim.prefix))
SH = load.bigmatrix(sprintf("%s.H", sim.prefix))
SP = load.bigmatrix(sprintf("%s.P", sim.prefix))

setwd(curdir)

err.prefix = sub("^(.+)\\.RData$", "\\1", basename(err.history))
EG = load.bigmatrix(sprintf("%s.G", err.prefix))
EH = load.bigmatrix(sprintf("%s.H", err.prefix))
EP = load.bigmatrix(sprintf("%s.P", err.prefix))



load(compare.file)



### genotype data

SG = as.matrix(SG)
EG = as.matrix(EG)

SAF = rowSums(SG) / (ncol(SG) * 2)
EAF = rowSums(EG) / (ncol(EG) * 2)


# focal site

d = rbind(data.table(type = "Missed rare variant",
										 index = diff.sim$index,
										 g0 = diff.sim$g0,
										 g1 = diff.sim$g1),
					data.table(type = "False rare variant",
										 index = diff.err$index,
										 g0 = diff.err$g0,
										 g1 = diff.err$g1),
					data.table(type = "True rare variant",
										 index = inter.err$index,
										 g0 = inter.err$g0,
										 g1 = inter.err$g1))

d$g0.in.err = EG[ cbind(d$index, d$g0) ]
d$g1.in.err = EG[ cbind(d$index, d$g1) ]

d$g0.in.sim = SG[ cbind(d$index, d$g0) ]
d$g1.in.sim = SG[ cbind(d$index, d$g1) ]

d$gg.sim = sprintf("%d%d", d$g0.in.sim, d$g1.in.sim)
d$gg.err = sprintf("%d%d", d$g0.in.err, d$g1.in.err)

i = which(d$gg.sim == "10"); if (length(i) > 0) d$gg.sim[i] = "01"
i = which(d$gg.sim == "20"); if (length(i) > 0) d$gg.sim[i] = "02"
i = which(d$gg.sim == "21"); if (length(i) > 0) d$gg.sim[i] = "12"

i = which(d$gg.err == "10"); if (length(i) > 0) d$gg.err[i] = "01"
i = which(d$gg.err == "20"); if (length(i) > 0) d$gg.err[i] = "02"
i = which(d$gg.err == "21"); if (length(i) > 0) d$gg.err[i] = "12"

d$from.to = sprintf("%s to %s", d$gg.sim, d$gg.err)


gg = ggplot(data = d) + 
	facet_grid(~type) +
	geom_bar(aes(from.to), colour="black") +
	scale_y_continuous(breaks = seq(0, 10e6, by=1e6*0.5), labels = seq(0, 10, by=0.5)) +
	theme_economist_white() +
	theme(legend.title = element_blank()) +
	xlab("Change in genotype pair at focal site") +
	ylab("Count (millions)") +
	ggtitle("Genotype pairs at focal rare variants")

ggsave(filename = "_plot.sim-compare.genotype.focal.pdf", plot = gg, width = 25, height = 10)



# breakpoints

d = rbind(data.table(type = "Detected on false rare variant",
										 g0 = c(diff.err$g0, diff.err$g0),
										 g1 = c(diff.err$g1, diff.err$g1),
										 brk = c(match(diff.err$g.brk.lhs, POS), match(diff.err$g.brk.rhs, POS)),
										 side = c(rep("LHS", nrow(diff.err)), rep("RHS", nrow(diff.err)))),
					data.table(type = "Detected on true rare variant",
										 g0 = c(inter.err$g0, inter.err$g0), 
										 g1 = c(inter.err$g1, inter.err$g1), 
										 brk = c(match(inter.err$g.brk.lhs, POS), match(inter.err$g.brk.rhs, POS)),
										 side = c(rep("LHS", nrow(inter.err)), rep("RHS", nrow(inter.err)))))

d$saf = SAF[d$brk]
d$eaf = EAF[d$brk]

d$g0.in.err = EG[ cbind(d$brk, d$g0) ]
d$g1.in.err = EG[ cbind(d$brk, d$g1) ]

d$g0.in.sim = SG[ cbind(d$brk, d$g0) ]
d$g1.in.sim = SG[ cbind(d$brk, d$g1) ]

d$gg.sim = sprintf("%d%d", d$g0.in.sim, d$g1.in.sim)
d$gg.err = sprintf("%d%d", d$g0.in.err, d$g1.in.err)

i = which(d$gg.sim == "10"); if (length(i) > 0) d$gg.sim[i] = "01"
i = which(d$gg.sim == "20"); if (length(i) > 0) d$gg.sim[i] = "02"
i = which(d$gg.sim == "21"); if (length(i) > 0) d$gg.sim[i] = "12"

i = which(d$gg.err == "10"); if (length(i) > 0) d$gg.err[i] = "01"
i = which(d$gg.err == "20"); if (length(i) > 0) d$gg.err[i] = "02"
i = which(d$gg.err == "21"); if (length(i) > 0) d$gg.err[i] = "12"

d$from.to = sprintf("%s to %s", d$gg.sim, d$gg.err)


gg = ggplot(data = d) + 
	facet_grid(~type) +
	geom_bar(aes(from.to, fill=side), colour="black", position="dodge") +
	scale_y_continuous(breaks = seq(0, 10e6, by=5e5), labels = seq(0, 10, by=0.5)) +
	scale_fill_economist() +
	theme_economist_white() +
	theme(legend.title = element_blank()) +
	xlab("Change in genotype pair at breakpoint site") +
	ylab("Count (millions)")+
	ggtitle("Genotype pairs at breakpoint sites")

ggsave(filename = "_plot.sim-compare.genotype.break-sides.pdf", plot = gg, width = 15, height = 10)


gg = ggplot(data = d) + 
	facet_grid(~type) +
	geom_bar(aes(from.to), colour="black", position="dodge") +
	scale_y_continuous(breaks = seq(0, 10e6, by=5e5), labels = seq(0, 10, by=0.5)) +
	theme_economist_white() +
	theme(legend.title = element_blank()) +
	xlab("Change in genotype pair at breakpoint site") +
	ylab("Count (millions)")+
	ggtitle("Genotype pairs at breakpoint sites")

ggsave(filename = "_plot.sim-compare.genotype.break.pdf", plot = gg, width = 15, height = 10)



# breakpoints by frequency

p = data.frame(brk = c(d$brk, d$brk), 
							 idv = c(d$g0, d$g1),
							 sim = c(d$g0.in.sim, d$g1.in.sim),
							 err = c(d$g0.in.err, d$g1.in.err),
							 saf = c(d$saf, d$saf),
							 stringsAsFactors = F)

p = unique(p)

splt = split(p, p$sim)
splt = lapply(splt, function(x) {
	
	brk = c(-1, unique(round((((0:250)/250)^exp(1)) * 5000)))
	idx = cut(round(x$saf * 5000), brk, include.lowest = T)
	
	z = split(x, idx)
	z = lapply(z, function(y) {
		if (nrow(y) == 0) return(NULL)
		n0 = length(which(y$err == 0)) 
		n1 = length(which(y$err == 1))
		n2 = length(which(y$err == 2))
		data.table(frq = median(y$saf), n0, n1, n2, p0 = n0 / nrow(y), p1 = n1 / nrow(y), p2 = n2 / nrow(y))
	})
	z = rbindlist(z)

	gt = sprintf("Simulated genotype is %d", x$sim[1])
	rbind(data.table(true.gt = gt, call.gt = "Error genotype is 0", count = z$n0, prop = z$p0, freq = z$frq), 
				data.table(true.gt = gt, call.gt = "Error genotype is 1", count = z$n1, prop = z$p1, freq = z$frq), 
				data.table(true.gt = gt, call.gt = "Error genotype is 2", count = z$n2, prop = z$p2, freq = z$frq) )
})
p = rbindlist(splt)


gg = ggplot(data = p) + 
	facet_grid(call.gt~true.gt) +
	geom_point(aes(x = freq, y = prop), colour = "orangered", alpha = 2/4) +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Allele frequency at breakpoint site") +
	ylab("Relative error proportion") +
	ggtitle("Genotypes at unique breakpoint sites")

ggsave(filename = "_plot.sim-compare.genotype.error-prop.pdf", plot = gg, width = 12, height = 12)


gg = ggplot(data = p) + 
	facet_grid(call.gt~true.gt) +
	geom_line(aes(x = freq, y = count), colour = "darkblue") +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Allele frequency at breakpoint site") +
	ylab("Relative error count") +
	ggtitle("Genotypes at unique breakpoint sites")

ggsave(filename = "_plot.sim-compare.genotype.error-count.pdf", plot = gg, width = 12, height = 12)







