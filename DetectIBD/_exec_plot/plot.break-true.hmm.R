#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("data.%s.RData", prefix))

load(sprintf("result.match_hmm_trunc.%s.RData", prefix))


DIST = c(0, cumsum(RATE[-(length(RATE))] * diff(POS) * 1e-6)) ### genetic distance
names(DIST) = as.character(POS)


del = which(match$g.wall | match$h.wall | match$p.wall | match$true.wall) # | match$hmm.wall)
if (length(del) > 0) {
	match = match[-del, ]
}

d = rbind(data.table(side="LHS", fk = match$fk, type2 = "DHG", type = "Discordant homozygote genotypes",
										 det  = match$pos - match$g.lhs,  tru  = match$pos - POS[match$true.lhs.idx]),
					data.table(side="LHS", fk = match$fk, type2 = "FGT, true haplotypes", type = " Four-gamete test, known haplotypes ",  
										 det = match$pos - match$h.lhs,   tru = match$pos - POS[match$true.lhs.idx]),
					data.table(side="LHS", fk = match$fk, type2 = "FGT, phased haplotypes", type = " Four-gamete test, phased haplotypes ", 
										 det = match$pos - match$p.lhs,   tru = match$pos - POS[match$true.lhs.idx]),
					#data.table(side="LHS", fk = match$fk, type2 = "HMM", type = "Hidden Markov Model, genotype data",   
					#det = match$pos - match$hmm.lhs, tru = match$pos - POS[match$true.lhs.idx]),
					data.table(side="RHS", fk = match$fk, type2 = "DHG", type = "Discordant homozygote genotypes",      
										 det = match$g.rhs - match$pos,   tru = POS[match$true.rhs.idx] - match$pos),
					data.table(side="RHS", fk = match$fk, type2 = "FGT, true haplotypes", type = " Four-gamete test, known haplotypes ",  
										 det = match$h.rhs - match$pos,   tru = POS[match$true.rhs.idx] - match$pos),
					data.table(side="RHS", fk = match$fk, type2 = "FGT, phased haplotypes", type = " Four-gamete test, phased haplotypes ", 
										 det = match$p.rhs - match$pos,   tru = POS[match$true.rhs.idx] - match$pos)
					#data.table(side="RHS", fk = match$fk, type2 = "HMM", type = "Hidden Markov Model, genotype data",   
					#det = match$hmm.rhs - match$pos, tru = POS[match$true.rhs.idx] - match$pos)
					)

del = which(d$tru == 0)
if (length(del) > 0)
	d = d[-del, ]


### RMSE
rmse = (d$det / 1e6) - (d$tru / 1e6)
rmse = split(rmse, list(d$fk, d$type))
rmse = sapply(rmse, function(r) {
	sqrt(mean(r^2))
})
rmse = data.table(tag = "- genotype error",
							 rmse = rmse,
							 fk = as.numeric(sub("^([0-9]+)\\.(.+)$", "\\1", names(rmse))),
							 type = sub("^([0-9]+)\\.(.+)$", "\\2", names(rmse)))

tmp = rmse

gg = ggplot(rmse) + 
	facet_grid(tag~.) + 
	geom_vline(xintercept = sort(unique(rmse$fk)), colour = "grey80") +
	geom_point(aes(x = fk, y = rmse, colour = type, shape = type), size = 2) +
	scale_x_continuous(breaks = sort(unique(rmse$fk)), labels = sort(unique(rmse$fk))) +
	theme_few() +
	theme(aspect.ratio = 1/3,
				legend.position = "top",
				legend.title = element_blank()) +
	ylab("RMSE (Mb)")

ggsave(gg, filename = "_plot.hmm.rmse_compare.pdf", width = 12, height = 8)



d$det = d$det + 1
d$tru = d$tru + 1



### hist, one-sided

d$map = (d$det / d$tru)

p = lapply(split(d, d$type), function(x) {
	z = density(x$map, bw = 0.05, n = 64, from=0, to=2.1)
	z = data.table(fk=x$fk[1], type=x$type[1], type2=x$type2[1], x = z$x, y = z$y / sum(z$y))
	z
})
p = rbindlist(p)

hist = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	#geom_point(aes(x=x, y=y), alpha=0.75, size=1) +
	geom_bar(aes(x=x, y=y), stat = "identity", alpha = 0.75, colour = "white", size = 0.25, fill = "black") +
	geom_vline(xintercept=c(1), colour="white", alpha=0.75) +
	geom_vline(xintercept=c(1), colour="grey30", linetype = "22") +
	#geom_hline(yintercept = 0) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels =  c(0, 0.5, 1, 1.5, 2)) +
	#scale_fill_manual(values = c(DHG="limegreen", FGT="sienna2", HMM="royalblue2")) +
	coord_cartesian(ylim = c(0, 0.19), expand = F) +
	theme_few() +
	theme(aspect.ratio=1,
				legend.title=element_blank(),
				legend.position = "top",
				strip.text = element_text(face = "bold"),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance to true breakpoint") +
	ylab("Density")

ggsave(filename = sprintf("__plot.break-true.hist.%s.pdf", prefix), plot = hist, width = 6, height = 10)
ggsave(filename = sprintf("__plot.break-true.hist.%s.png", prefix), plot = hist, width = 6, height = 10)



### mapped

d$map = (d$det / d$tru)
z = which(d$side == "LHS")
d$map[z] = d$map[z] * -1


p = lapply(split(d, d$type), function(x) {
	lhs = density(x$map, bw = 0.05, n = 81, from=-2.1, to=0)
	lhs = data.table(fk=x$fk[1], type=x$type[1], type2=x$type2[1], x = lhs$x, y = lhs$y / sum(lhs$y))
	rhs = density(x$map, bw = 0.05, n = 81, from=0, to=2.1)
	rhs = data.table(fk=x$fk[1], type=x$type[1], type2=x$type2[1], x = rhs$x, y = rhs$y / sum(rhs$y))
	rbind(lhs, rhs)
})
p = rbindlist(p)


mapped = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	geom_bar(aes(x=x, y=y), stat = "identity", alpha = 0.75, colour = "white", size = 0.25, fill = "black") +
	geom_vline(xintercept=c(-1, 0, 1), colour="white", alpha=0.75) +
	geom_vline(xintercept=c(-1, 0, 1), colour="grey30", linetype = "22") +
	#geom_hline(yintercept = 0) +
	scale_x_continuous(breaks = c(2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2)) +
	coord_cartesian(xlim=c(-2.1, 2.1), ylim = c(0, 0.185), expand = F) +
	# geom_vline(xintercept=c(-1, 0, 1), colour="grey50", linetype = "22") +
	# geom_line(aes(x=x, y=y), alpha=0.75, size=1) +
	annotate("text", label = "LHS", x = -0.25, y = 0.175, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x =  0.25, y = 0.175, size = 3.5, colour="grey20") +
	# scale_x_continuous(breaks = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3), labels = c(3, 2, 1, 0.5, 0, 0.5, 1, 2, 3)) +
	theme_few() +
	theme(aspect.ratio=1/2,
				legend.title=element_blank(),
				legend.position = "top",
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance between focal site and detected breakpoint") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true.hmm.mapped.pdf", prefix), plot = mapped, width = 6, height = 10)




### heatmap, one-sided

brk = seq(log10(1), log10(100e6), length.out = 201)

p = lapply(split(d, d$type), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


p$Var1 = p$Var1 + 1


tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


heat = ggplot(data = p) +
	facet_wrap(~type, ncol = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(1), intercept = c(-1), colour="black", alpha=0.25) +
	#scale_fill_gradient2(low = "white", mid = "gold", high = "darkblue", na.value = "grey75", midpoint = max(log10(p$x)) * 0.25, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	#scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey75", breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_fill_gradient2(low = "grey70", mid = "beige", high = "darkgreen", na.value = "grey60", midpoint = max(log10(p$x)) * 0.5, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1/1,
				legend.title = element_blank(),
				axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")

ggsave(filename = sprintf("_plot.break-true.hmm.heat.pdf", prefix), plot = heat, width = 12, height = 12)
ggsave(filename = sprintf("_plot.break-true.hmm.heat.png", prefix), plot = heat, width = 12, height = 12)






### scatter

brk = seq(log10(1), log10(100e6), length.out = 201)

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
				strip.text = element_text(face = "bold"),
				legend.title = element_blank()) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")

ggsave(filename = sprintf("__plot.break-true.scatter.%s.pdf", prefix), plot = scatter, width = 10, height = 10)
ggsave(filename = sprintf("__plot.break-true.scatter.%s.png", prefix), plot = scatter, width = 10, height = 10)





###
### length
###


d = rbind(data.table(fk = match$fk, type = "Discordant homozygote genotypes",  det = (match$g.rhs - match$g.lhs),     tru = (match$true.rhs.pos - match$true.lhs.pos)),
					data.table(fk = match$fk, type = " Four-gamete test, known haplotypes ",  det = (match$h.rhs - match$h.lhs),     tru = (match$true.rhs.pos - match$true.lhs.pos)),
					data.table(fk = match$fk, type = " Four-gamete test, phased haplotypes ", det = (match$p.rhs - match$p.lhs),     tru = (match$true.rhs.pos - match$true.lhs.pos))
					#data.table(fk = match$fk, type = "Hidden Markov Model, genotype data",   det = (match$hmm.rhs - match$hmm.lhs), tru = (match$true.rhs.pos - match$true.lhs.pos)))
)

del = which(d$tru == 0)
if (length(del) > 0)
	d = d[-del, ]

#d$det = d$det + 1
#d$tru = d$tru + 1


# mapped

d$map = (d$det / d$tru)


p = lapply(split(d, d$type), function(x) {
	tmp = density(x$map, bw = 0.005, n = 1024*10, from=0, to=2)
	data.table(fk=x$fk[1], type=x$type[1], x = tmp$x, y = tmp$y / sum(tmp$y))
})
p = rbindlist(p)


mapped = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	geom_vline(xintercept=c(1), colour="grey50", linetype = "22") +
	geom_line(aes(x=x, y=y), alpha=0.75, size=1) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = c(0, 0.5, 1, 1.5, 2)) +
	theme_few() +
	theme(aspect.ratio=1/2,
				legend.title=element_blank(),
				legend.position = "top",
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Detected segment length relative to true length") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true.hmm.length-mapped.pdf", prefix), plot = mapped, width = 10, height = 14)


### hist

d$map = (d$det / d$tru)

p = lapply(split(d, d$type), function(x) {
	z = density(x$map, bw = 0.05, n = 64, from=0, to=2.1)
	z = data.table(fk=x$fk[1], type=x$type[1], x = z$x, y = z$y / sum(z$y))
	z
})
p = rbindlist(p)

hist = ggplot(data=p) +
	facet_wrap(~type, nrow = 1) +
	#geom_point(aes(x=x, y=y), alpha=0.75, size=1) +
	geom_bar(aes(x=x, y=y), stat = "identity", alpha = 0.75, colour = "white", size = 0.25, fill = "black") +
	geom_vline(xintercept=c(1), colour="white", alpha=0.75) +
	geom_vline(xintercept=c(1), colour="grey30", linetype = "22") +
	geom_hline(yintercept = 0) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels =  c(0, 0.5, 1, 1.5, 2)) +
	#scale_fill_manual(values = c(DHG="limegreen", FGT="sienna2", HMM="royalblue2")) +
	coord_cartesian(ylim = c(0, 0.185)) +
	theme_few() +
	theme(aspect.ratio=1,
				legend.title=element_blank(),
				legend.position = "top",
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Detected segment length relative to true length") +
	ylab("Density")

ggsave(filename = sprintf("_plot.break-true.hmm.hist-length.pdf", prefix), plot = hist, width = 12, height = 12)
ggsave(filename = sprintf("_plot.break-true.hmm.hist-length.png", prefix), plot = hist, width = 12, height = 12)





# scatter

brk = seq(log10(1), log10(100e6), length.out = 201)

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


scatter = ggplot(data = p) +
	facet_wrap(~type, ncol = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = 1, intercept = 0, colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = max(log10(p$x)) * 0.6, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold"),
				legend.title = element_blank()) +
	ylab("Detected segment length") +
	xlab("True segment length")

ggsave(filename = sprintf("__plot.break-true.length-scatter.%s.pdf", prefix), plot = scatter, width = 6, height = 10)
ggsave(filename = sprintf("__plot.break-true.length-scatter.%s.png", prefix), plot = scatter, width = 6, height = 10)



### heatmap

brk = seq(log10(1), log10(100e6), length.out = 201)

p = lapply(split(d, d$type), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


p$Var1 = p$Var1 + 1


tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


heat = ggplot(data = p) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = c(1), intercept = c(-1), colour="black", alpha=0.25) +
	#scale_fill_gradient2(low = "white", mid = "gold", high = "darkblue", na.value = "grey75", midpoint = max(log10(p$x)) * 0.25, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	#scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey75", breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_fill_gradient2(low = "grey70", mid = "beige", high = "darkgreen", na.value = "grey60", midpoint = max(log10(p$x)) * 0.5, breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1/1,
				legend.title = element_blank(),
				axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
	ylab("Detected segment length") +
	xlab("True segment length")

ggsave(filename = sprintf("_plot.break-true.hmm.heat-length.pdf", prefix), plot = heat, width = 12, height = 12)
ggsave(filename = sprintf("_plot.break-true.hmm.heat-length.png", prefix), plot = heat, width = 12, height = 12)











