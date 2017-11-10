#
# plot error profile
#

library(ggplot2)
library(grid)


pro.file = "profile.ALL.RData"


load(pro.file)


del = which(p$maf == 0 & !p$error)
if (length(del) > 0) {
	p = p[-del, ]
}



col = c("Truth 0, called 0" = "peachpuff",
				"Truth 0, called 1" = "orange",
				"Truth 0, called 2" = "orangered",
				"Truth 1, called 0" = "green2",
				"Truth 1, called 1" = "darkseagreen1",
				"Truth 1, called 2" = "seagreen",
				"Truth 2, called 0" = "dodgerblue4",
				"Truth 2, called 1" = "deepskyblue",
				"Truth 2, called 2" = "lightblue1")



d = p

brk = seq(0, 0.5, length.out = 101) # c(0, exp(seq(log(0.0005), log(0.5), length.out = 100))) # 

d$bin = cut(d$maf, breaks = brk, include.lowest = T)

splt = split(d, d$bin)
splt = lapply(splt, function(x) {
	sub = split(x, x$type)
	sub = lapply(sub, function(y, n) {
		z = y[1, ]
		z$MAF = mean(y$maf)
		z$count = nrow(y)
		z$percent = nrow(y) / n * 100
		z
	}, nrow(x))
	sub = Reduce(rbind, sub)
	if (nrow(x) == 0) return(NULL)
	sub$nbin = nrow(x)
	sub
})
d = Reduce(rbind, splt)



d.prop = d
ymax = 100
### OR ###
d.prop = d[d$error,]
ymax = 2.5

gg = ggplot(data = d.prop) + 
	geom_bar(aes(x=as.numeric(bin)-0.5, y=percent, fill=label), stat="identity") + 
	scale_fill_manual(values = col) +
	scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, ymax, length.out = 11)) +
	coord_cartesian(ylim = c(0, ymax)) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank()) +
	xlab("MAF bin (%)") +
	ylab("Relative Proportion (%)")

ggsave(sprintf("_plot.prop_complete_together_stretched.%s.pdf", pro.file), gg, width = 12, height = 24)
ggsave(sprintf("_plot.prop_complete_together_normal.%s.pdf", pro.file), gg, width = 12, height = 12)
### OR ###
ggsave(sprintf("_plot.prop_error_together.%s.pdf", pro.file), gg, width = 12, height = 12)



d.prop = d
ymax = 100
### OR ###
d.prop = d[d$error,]
ymax = 1.5

gg = ggplot(data = d.prop) + 
	geom_line(aes(x=as.numeric(bin)-0.5, y=percent, colour=label), size=1) + 
	scale_colour_manual(values = col) +
	scale_x_continuous(expand = c(0.005, 0.005), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_continuous(expand = c(0.005, 0.005), breaks = seq(0, ymax, length.out = 7)) +
	coord_cartesian(ylim = c(0, ymax)) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank()) +
	xlab("MAF bin (%)") +
	ylab("Relative Proportion (%)")

ggsave(sprintf("_plot.lines_prop_complete_together_stretched.%s.pdf", pro.file), gg, width = 12, height = 24)
ggsave(sprintf("_plot.lines_prop_complete_together_normal.%s.pdf", pro.file), gg, width = 12, height = 12)
### OR ###
ggsave(sprintf("_plot.lines_prop_error_together.%s.pdf", pro.file), gg, width = 12, height = 12)



d.size = data.frame(bin = d$bin, nbin = d$nbin / 1e03)
d.size = unique(d.size)
if (any(d.size$nbin < 1)) {
	d.size$nbin[which(d.size$nbin < 1)] = 1
}

gg = ggplot(data = d.size) + 
	geom_bar(aes(x=as.numeric(bin)-0.5, y=nbin), alpha=0.5, stat="identity") + 
	scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_log10(expand = c(0.0025, 0.0025), breaks = 10^(0:8) ) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank()) +
	xlab("MAF bin (%)") +
	ylab("Number of markers (thousands), log-scale")

ggsave(sprintf("_plot.size_complete_together.%s.pdf", pro.file), gg, width = 12, height = 6)





d = p

brk = seq(0, 0.5, length.out = 101) # c(0, exp(seq(log(0.0005), log(0.5), length.out = 100))) # 

d$bin = cut(d$maf, breaks = brk, include.lowest = T)

splt = split(d, list(d$bin, d$ref.gt))
splt = lapply(splt, function(x) {
	sub = split(x, x$type)
	sub = lapply(sub, function(y, n) {
		z = y[1, ]
		z$MAF = mean(y$maf)
		z$count = nrow(y)
		z$percent = nrow(y) / n * 100
		z
	}, nrow(x))
	sub = Reduce(rbind, sub)
	if (nrow(x) == 0) return(NULL)
	sub$nbin = nrow(x)
	sub
})
d = Reduce(rbind, splt)


d$true.gt = NA
d$true.gt[ which(d$ref.gt == 0) ] = "True genotype = 0"
d$true.gt[ which(d$ref.gt == 1) ] = "True genotype = 1"
d$true.gt[ which(d$ref.gt == 2) ] = "True genotype = 2"


d.prop = d
ymax = 100
### OR ###
d.prop = d[d$error,]
ymax = 8

gg = ggplot(data = d.prop) + 
	facet_wrap(~true.gt) +
	geom_bar(aes(x=as.numeric(bin)-0.5, y=percent, fill=label), stat="identity") + 
	scale_fill_manual(values = col) +
	scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, ymax, length.out = 21)) +
	coord_cartesian(ylim = c(0, ymax)) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank(),
				strip.background = element_rect(fill = NA, colour = NA),
				panel.border = element_rect(fill = NA, colour = "black"),
				panel.margin = unit(0.25, units = "inches")) +
	xlab("MAF bin (%)") +
	ylab("Relative Proportion (%)")

ggsave(sprintf("_plot.prop_complete_bytruegt_stretched.%s.pdf", pro.file), gg, width = 12, height = 24)
ggsave(sprintf("_plot.prop_complete_bytruegt_normal.%s.pdf", pro.file), gg, width = 12, height = 12)
### OR ###
ggsave(sprintf("_plot.prop_error_bytruegt.%s.pdf", pro.file), gg, width = 12, height = 12)



d.prop = d
ymax = 100
### OR ###
d.prop = d[d$error,]
ymax = 8

gg = ggplot(data = d.prop) + 
	facet_wrap(~true.gt) +
	geom_line(aes(x=as.numeric(bin)-0.5, y=percent, colour=label), size=1) + 
	scale_colour_manual(values = col) +
	scale_x_continuous(expand = c(0.005, 0.005), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_continuous(expand = c(0.005, 0.005), breaks = seq(0, ymax, length.out = 21)) +
	coord_cartesian(ylim = c(0, ymax)) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank(),
				strip.background = element_rect(fill = NA, colour = NA),
				panel.border = element_rect(fill = NA, colour = "black"),
				panel.margin = unit(0.25, units = "inches")) +
	xlab("MAF bin (%)") +
	ylab("Relative Proportion (%)")

ggsave(sprintf("_plot.lines_prop_complete_bytruegt_stretched.%s.pdf", pro.file), gg, width = 12, height = 24)
ggsave(sprintf("_plot.lines_prop_complete_bytruegt_normal.%s.pdf", pro.file), gg, width = 12, height = 12)
### OR ###
ggsave(sprintf("_plot.lines_prop_error_bytruegt.%s.pdf", pro.file), gg, width = 12, height = 12)



d.size = data.frame(bin = d$bin, nbin = d$nbin / 1e03, true.gt = d$true.gt)
d.size = unique(d.size)
if (any(d.size$nbin < 1)) {
	d.size$nbin[which(d.size$nbin < 1)] = 1
}

gg = ggplot(data = d.size) + 
	facet_wrap(~true.gt) +
	geom_bar(aes(x=as.numeric(bin)-0.5, y=nbin), alpha=0.5, stat="identity") + 
	scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
	scale_y_log10(expand = c(0.0025, 0.0025), breaks = 10^(0:8) ) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				panel.grid.major.x = element_blank(),
				legend.title = element_blank(),
				strip.background = element_rect(fill = NA, colour = NA),
				panel.border = element_rect(fill = NA, colour = "black"),
				panel.margin = unit(0.25, units = "inches")) +
	xlab("MAF bin (%)") +
	ylab("Number of markers (thousands), log-scale")

ggsave(sprintf("_plot.size_complete_bytruegt.%s.pdf", pro.file), gg, width = 12, height = 6)





gt = as.character(0:2)
m = matrix(0, nrow = 3, ncol = 3, dimnames = list(gt, gt))
for (tgt in gt) {
	for (cgt in gt) {
		m[tgt, cgt] = length(which(p$ref.gt == tgt & p$alt.gt == cgt))
	}
}
m = m / apply(m, 1, sum)















