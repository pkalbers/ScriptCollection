

library(data.table)
library(ggplot2)
library(ggthemes)


load("../frq.RData")

frq = split(frq, frq$chr)
num = lapply(frq, function(x) {
	z = table(x$alt.frq)
	data.table(f = as.numeric(names(z)), n = as.vector(z))
})


save(num, file = "frq.table.R")


f = sort(unique(unlist(sapply(num, function(x) x$f))))


d = NULL

for (i in 1:22) {
	tmp = data.table(f = f, n = 0, chr = sprintf("% 3d", i))
	sub = num[[i]]
	x = match(sub$f, as.character(f))
	tmp$n[x] = sub$n
	d = rbind(d, tmp)
}


del = which(d$f < 0.0002 | d$n == 0)
d = d[-del,]


d = split(d, d$chr)
d = lapply(d, function(d) {
	d = d[order(d$f),]

	a = diff(log10(c(d$n[1], d$n)))
	b = diff(log10(c(d$n, d$n[nrow(d)])))
	z = abs(a - b)
	del = which(d$f > 0.005 & z < quantile(z, 0.95))
	del = del[which(1:length(del) %% 5 != 0)]

	if (length(del) > 0) d = d[-del, ]

	d
})
d = rbindlist(d)


d$chr = factor(d$chr, levels = unique(d$chr), ordered = T)


xtik = c(0.0005, 0.001,
				 seq(0.002, 0.01, length.out = 5),
				 seq(0.02, 0.1, length.out = 5),
				 seq(0.2, 1, length.out = 5))
#xtik = c(c(0.0005,0.005,0.05,0.5), c(0.001,0.01,0.1,1))
xmin = unique(c(seq(0.0001, 0.001, length.out = 10),
								seq(0.001, 0.01, length.out = 10),
								seq(0.01, 0.1, length.out = 10),
								seq(0.1, 1, length.out = 10),
								seq(1, 10, length.out = 10)))
ymin = c(xmin * 1e4, xmin * 1e6)


gg = ggplot(d) +
	geom_line(aes(x = f, y = n, color = chr)) +
	scale_x_log10(breaks = xtik, labels = trimws(format(xtik * 100, scientific = F, drop0trailing = T)), minor_breaks = xmin) +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), minor_breaks = ymin,
								labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000")) +
	scale_color_manual(values = colorRampPalette(rev(c("red4", "orange", "yellow", "darkgreen", "cyan", "royalblue", "navy")))(22), name = "Chromosome") +
	coord_cartesian(expand = F, ylim = c(0.975, 1.5e6), xlim = c(0.000395, 1.015)) +
	#theme_classic() +
	theme(aspect.ratio = 9/16,
				panel.background = element_rect(fill = "grey85"),
				legend.key.height = unit(0.4, "cm"),
				legend.position = c(1, 1),
				panel.grid.major = element_line(size = 1/3),
				panel.grid.minor = element_line(size = 1/3),
				legend.justification = c(1.025, 1.01),
				panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
	xlab("Allele frequency (%)") + ylab("Count") +
	guides(colour = guide_legend(override.aes = list(size=2)))
gg

ggsave(gg, filename = "allelefreq_1kg.pdf", width = 9, height = 5)




