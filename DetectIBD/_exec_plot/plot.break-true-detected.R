#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("data.%s.RData", prefix))

load(sprintf("result.match.%s.RData", prefix))
### OR:
load(sprintf("result.beagle.%s.RData", prefix))


del = which(match$g.wall | match$h.wall | match$p.wall | match$true.wall)
if (length(del) > 0) {
	match = match[-del, ]
}

d = rbind(data.table(side="LHS", fk = match$fk, type = "Discordant Homozygote Genotype Test",  det = match$pos - match$g.lhs, tru = match$pos - POS[match$true.lhs.idx]),
					data.table(side="LHS", fk = match$fk, type = "Four-Gametes Test, known haplotypes",  det = match$pos - match$h.lhs, tru = match$pos - POS[match$true.lhs.idx]),
					data.table(side="LHS", fk = match$fk, type = "Four-Gametes Test, phased haplotypes", det = match$pos - match$p.lhs, tru = match$pos - POS[match$true.lhs.idx]),
					data.table(side="RHS", fk = match$fk, type = "Discordant Homozygote Genotype Test",  det = match$g.rhs - match$pos, tru = POS[match$true.rhs.idx] - match$pos),
					data.table(side="RHS", fk = match$fk, type = "Four-Gametes Test, known haplotypes",  det = match$h.rhs - match$pos, tru = POS[match$true.rhs.idx] - match$pos),
					data.table(side="RHS", fk = match$fk, type = "Four-Gametes Test, phased haplotypes", det = match$p.rhs - match$pos, tru = POS[match$true.rhs.idx] - match$pos))

### OR:

del = which(match.H$wall)
if (length(del) > 0) {
	match.H = match.H[-del, ]
}
del = which(match.P$wall)
if (length(del) > 0) {
	match.P = match.P[-del, ]
}

d = rbind(data.table(side="LHS", fk = match.H$fk, type = "Beagle IBD, known haplotypes",   det = match.H$position - match.H$beagle.lhs.position, tru = match.H$position - POS[match.H$lhs.index]),
					data.table(side="LHS", fk = match.P$fk, type = "Beagle IBD, phased haplotypes",  det = match.P$position - match.P$beagle.lhs.position, tru = match.P$position - POS[match.P$lhs.index]),
					data.table(side="RHS", fk = match.H$fk, type = "Beagle IBD, known haplotypes",   det = match.H$beagle.rhs.position - match.H$position, tru = POS[match.H$rhs.index] - match.H$position),
					data.table(side="RHS", fk = match.P$fk, type = "Beagle IBD, phased haplotypes",  det = match.P$beagle.rhs.position - match.P$position, tru = POS[match.P$rhs.index] - match.P$position))



d$det = d$det + 1
d$tru = d$tru + 1


### mapped

d$map = (d$det / d$tru)
z = which(d$side == "LHS")
d$map[z] = d$map[z] * -1


p = lapply(split(d, d$type), function(x) {
	lhs = density(x$map, bw = 0.005, n = 1024*10, from=-3.5, to=0)
	lhs = data.table(fk=x$fk[1], type=x$type[1], x = lhs$x, y = lhs$y / sum(lhs$y))
	rhs = density(x$map, bw = 0.005, n = 1024*10, from=0, to=3.5)
	rhs = data.table(fk=x$fk[1], type=x$type[1], x = rhs$x, y = rhs$y / sum(rhs$y))
	rbind(lhs, rhs)
})
p = rbindlist(p)


mapped = ggplot(data=p) + 
	facet_wrap(~type, ncol = 1) + 
	geom_vline(xintercept=c(-1, 0, 1), colour="grey50", linetype = "22") +
	geom_line(aes(x=x, y=y), alpha=0.75, size=1) + 
	annotate("text", label = "LHS", x = -3.5, y = max(p$y) * 0.975, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x =  3.5, y = max(p$y) * 0.975, size = 3.5, colour="grey20") +
	scale_x_continuous(breaks = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3), labels = c(3, 2, 1, 0.5, 0, 0.5, 1, 2, 3)) +
	theme_few() +
	theme(aspect.ratio=1/2,
				legend.title=element_blank(),
				legend.position = "top", 
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Detected breakpoint distance to focal site relative to true breakpoint") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true-detected.naive.mapped.pdf", prefix), plot = mapped, width = 8, height = 12)
### OR:
ggsave(filename = sprintf("_plot.break-true-detected.beagle.mapped.pdf", prefix), plot = mapped, width = 8, height = 8)



### scatter

brk = seq(log10(1), log10(100e6), length.out = 251)

p = lapply(split(d, list(d$type, d$side)), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	#mat$fk  = x$fk[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


p$Var1 = p$Var1 + 1
x = (p$side == "LHS")
p$Var1[x] = p$Var1[x] * -1


tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T)) 


scatter = ggplot(data = p) + 
	facet_wrap(~type, ncol = 1) + 
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(-1, 1), intercept = c(-1, -1), colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1/2,
				legend.title = element_blank()) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")

ggsave(filename = sprintf("_plot.break-true-detected.naive.scatter.pdf", prefix), plot = scatter, width = 8, height = 12)
### OR:
ggsave(filename = sprintf("_plot.break-true-detected.beagle.scatter.pdf", prefix), plot = scatter, width = 8, height = 8)



### length

d = rbind(data.table(fk = match$fk, type = "Discordant Homozygote Genotype Test",  det = (match$g.rhs - match$g.lhs), tru = (POS[match$true.rhs.idx] - POS[match$true.lhs.idx])),
					data.table(fk = match$fk, type = "Four-Gametes Test, known haplotypes",  det = (match$h.rhs - match$h.lhs), tru = (POS[match$true.rhs.idx] - POS[match$true.lhs.idx])),
					data.table(fk = match$fk, type = "Four-Gametes Test, phased haplotypes", det = (match$p.rhs - match$p.lhs), tru = (POS[match$true.rhs.idx] - POS[match$true.lhs.idx])))

### OR:

d = rbind(data.table(fk = match.H$fk, type = "Beagle IBD, known haplotypes",   det = (match.H$beagle.rhs.position - match.H$beagle.lhs.position), tru = (POS[match.H$rhs.index] - POS[match.H$lhs.index])),
					data.table(fk = match.P$fk, type = "Beagle IBD, phased haplotypes",  det = (match.P$beagle.rhs.position - match.P$beagle.lhs.position), tru = (POS[match.P$rhs.index] - POS[match.P$lhs.index])))



brk = seq(log10(1), log10(100e6), length.out = 251)

p = lapply(split(d, d$type), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)



tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T)) 


len = ggplot(data = p) + 
	facet_wrap(~type, ncol = 1) + 
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = 1, intercept = 0, colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	ylab("Detected segment length") +
	xlab("True segment length")

ggsave(filename = sprintf("_plot.break-true-detected.naive.length.pdf", prefix), plot = len, width = 8, height = 12)
### OR:
ggsave(filename = sprintf("_plot.break-true-detected.beagle.length.pdf", prefix), plot = len, width = 8, height = 8)






d = data.table(fk = truth$fk, time = truth$time, size = truth$rhs.position - truth$lhs.position)

ggplot(d[sample(nrow(d), 1e5),]) + geom_point(aes(x = size, y = time, colour = fk), alpha = 0.1) + scale_x_log10() + scale_y_log10()





