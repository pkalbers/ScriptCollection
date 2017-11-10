#
# plot detected vs truth
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)
library(cowplot)


args = commandArgs(T)

prefix = args[1]             # file prefix


load(sprintf("result_hmm.match.%s.RData", prefix))

hmm = match

load(sprintf("result.match.%s.RData", prefix))


key.hmm = sprintf("%d %d %d", hmm$index, hmm$h0, hmm$h1)
key.mat = sprintf("%d %d %d", match$index, match$h0, match$h1)
key = intersect(key.hmm, key.mat)

i.hmm = match(key, key.hmm)
i.mat = match(key, key.mat)

hmm = hmm[i.hmm]
match = match[i.mat]

match = cbind(match, as.data.frame(hmm)[ grep("^m\\.+", names(hmm)) ])



d = rbind(data.table(side="LHS", fk = match$fk, type = "Discordant Homozygote Genotype Test",  det = match$g.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Discordant Homozygote Genotype Test",  det = match$g.brk.dist.rhs, tru = match$t.brk.dist.rhs),
					data.table(side="LHS", fk = match$fk, type = "Four-Gametes Test, known haplotypes",  det = match$h.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Four-Gametes Test, known haplotypes",  det = match$h.brk.dist.rhs, tru = match$t.brk.dist.rhs),
					data.table(side="LHS", fk = match$fk, type = "Four-Gametes Test, phased haplotypes", det = match$p.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Four-Gametes Test, phased haplotypes", det = match$p.brk.dist.rhs, tru = match$t.brk.dist.rhs),
					data.table(side="LHS", fk = match$fk, type = "Hidden Markov Model, genotype data",   det = match$m.brk.dist.lhs, tru = match$t.brk.dist.lhs),
					data.table(side="RHS", fk = match$fk, type = "Hidden Markov Model, genotype data",   det = match$m.brk.dist.rhs, tru = match$t.brk.dist.rhs))


### mapped

d$map = (d$det / d$tru)
z = which(d$side == "LHS")
d$map[z] = d$map[z] * -1


p = lapply(split(d, d$type), function(x) {
	lhs = density(x$map, bw = 0.005, n = 1024*10, from=-2.5, to=0)
	lhs = data.table(fk=x$fk[1], type=x$type[1], x = lhs$x, y = lhs$y / sum(lhs$y))
	rhs = density(x$map, bw = 0.005, n = 1024*10, from=0, to=2.5)
	rhs = data.table(fk=x$fk[1], type=x$type[1], x = rhs$x, y = rhs$y / sum(rhs$y))
	rbind(lhs, rhs)
})
p = rbindlist(p)


mapped = ggplot(data=p) + 
	facet_wrap(~type, ncol = 1) + 
	geom_vline(xintercept=c(-1, 0, 1), colour="grey50", linetype = "22") +
	geom_line(aes(x=x, y=y), alpha=0.75, size=1) + 
	scale_x_continuous(breaks = -5:5) +
	theme_few() +
	theme(aspect.ratio=1/2,
				legend.title=element_blank(),
				legend.position = "top", 
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Relative physical distance from focal site to true breakpoints") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true-detected.hmm.mapped.pdf", prefix), plot = mapped, width = 10, height = 5)



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
tik = as.numeric(cut(log10(c(100, 10000, 1e6, 100e6)), breaks = brk, include.lowest = T)) 


scatter = ggplot(data = p) + 
	facet_wrap(~type, ncol = 1) + 
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(-1, 1), intercept = c(-1, -1), colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1/2) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance (IBD)")

ggsave(filename = sprintf("_plot.break-true-detected.hmm.scatter.pdf", prefix), plot = scatter, width = 10, height = 5)




gg = plot_grid(mapped, scatter, ncol = 2, rel_widths = c(8.5/10, 1))

ggsave(filename = sprintf("_plot.break-true-detected.all.pdf", prefix), plot = gg, width = 12, height = 12)







