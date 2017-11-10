
library(data.table)
library(ggplot2)
library(ggthemes)


load("../data.OutOfAfricaHapMap20.RData")

load("phase.rawdata.1276.RData")



d = coord
d = split(d, d$other)
d = lapply(d, function(x) {
	if (nrow(x) == 1) return (x)
	key = sprintf("%d %d", x$true.beg, x$true.end)
	del = which(duplicated(key))
	if (length(del) > 0) {
		x = x[-del,]
	}
	x
})
d = rbindlist(d)


cov = list()

chr = rep(0, length(POS))

fks = as.character(sort(unique(d$fk)))

d = split(d, d$fk)

for (fk in fks) {
	print(fk)
	k = d[[fk]]
	for (i in 1:nrow(k)) {
		rng = (k$true.beg[i]):(k$true.end[i])
		chr[rng] = chr[rng] + 1
	}
	cov[[fk]] = chr
}

d = rbindlist(d)
p = lapply(cov, function(x) data.table(idx = 1:length(POS), pos = POS, cov = x))

for (fk in fks) {
	q = p[[fk]]
	q$fk = as.numeric(fk)
	x = which(diff(q$cov) != 0)
	x = sort(unique(c(1, x, x+1, nrow(q))))
	p[[fk]] = q[x, ]
}
p = rbindlist(p)




d = p
d$fk = factor(d$fk, levels = as.character(rev(2:25)), ordered = T)

gg = ggplot(d) + 
	geom_ribbon(aes(x=pos, ymax=cov, ymin=0, fill = fk)) +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6), expand = c(0,0)) +
	scale_y_continuous(breaks = seq(0, 500, by = 10), expand = c(0,0)) +
	scale_fill_manual(values = colorRampPalette(c(colorRampPalette(c("grey90", "purple"))(3)[2], 
																								colorRampPalette(c("grey70", "navy"))(3)[2], 
																								#colorRampPalette(c("grey50", "green"))(3)[2], 
																								colorRampPalette(c("grey20", "red"))(3)[2], 
																								colorRampPalette(c("grey95", "yellow"))(3)[2]))(24)) +
	coord_cartesian(ylim = c(1, 55)) +
	theme_few() + 
	theme(legend.title = element_blank(),
				panel.border = element_rect(fill = NA, colour = "black", size = 1),
				panel.grid.major.y = element_line(colour = "grey80", size = 0.5),
				legend.key.size = unit(0.5, "cm"),
				legend.direction = "horizontal",
				legend.position = "bottom") +
	xlab("Physical position (Mb)") +
	ylab("Cumulative coverage") + 
	guides(fill = guide_legend(nrow = 1, reverse = TRUE))

ggsave(gg, filename = "plot.fk_coverage_example.pdf", width = 10, height = 5)



