#
# plot detected vs truth
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


args = commandArgs(T)

prefix = args[1]             # file prefix


load(sprintf("result.match.%s.RData", prefix))



d = rbind(data.table(side="LHS", fk = match$fk, type = "Genotype",         det = match$g.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Genotype",         det = match$g.brk.dist.rhs, tru = match$t.brk.dist.rhs),
					data.table(side="LHS", fk = match$fk, type = "Haplotype",        det = match$h.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Haplotype",        det = match$h.brk.dist.rhs, tru = match$t.brk.dist.rhs),
					data.table(side="LHS", fk = match$fk, type = "Phased haplotype", det = match$p.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Phased haplotype", det = match$p.brk.dist.rhs, tru = match$t.brk.dist.rhs))



### mapped

d$map = (d$det / d$tru)
z = which(d$side == "LHS")
d$map[z] = d$map[z] * -1


p = lapply(split(d, d$type), function(x) {
	lhs = density(x$map, bw = 0.005, n = 1024*10, from=-5, to=0)
	lhs = data.table(type=x$type[1], x = lhs$x, y = lhs$y / sum(lhs$y))
	rhs = density(x$map, bw = 0.005, n = 1024*10, from=0, to=5)
	rhs = data.table(type=x$type[1], x = rhs$x, y = rhs$y / sum(rhs$y))
	rbind(lhs, rhs)
})
p = rbindlist(p)


gg = ggplot(data=p) + 
	geom_vline(xintercept=c(-1, 0, 1), colour="grey") +
	geom_line(aes(x=x, y=y, colour=type), alpha=0.75, size=1) + 
	scale_x_continuous(breaks = -5:5) +
	theme_few() +
	theme(aspect.ratio=1/3,
				legend.title=element_blank(),
				legend.position = "top") +
	xlab("Relative physical distance") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true-detected.mapped.pdf", prefix), plot = gg, width = 15, height = 5)


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


tag = c("100b", "10Kb", "1Mb", "100Mb")
tik = as.numeric(cut(log10(c(100, 10000, 1e6, 100e6)), breaks = brk, include.lowest = T)) + 1


gg = ggplot(data = p) + 
	facet_grid(.~type) + 
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

ggsave(filename = sprintf("_plot.break-true-detected.scatter.pdf", prefix), plot = gg, width = 15, height = 5)




