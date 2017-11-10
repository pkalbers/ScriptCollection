
library(data.table)
library(ggplot2)
library(ggthemes)


q = (0:500)/500
p = 1 - q

e = rbind(data.table(gt = " 0, 0", x = p, y = p^4, type = " Non-IBD "),
					data.table(gt = " 0, 1", x = p, y = 4*p^3*q, type = " Non-IBD "),
					data.table(gt = " 0, 2", x = p, y = 2*p^2*q^2, type = " Non-IBD "),
					data.table(gt = " 1, 1", x = p, y = 4*p^2*q^2, type = " Non-IBD "),
					data.table(gt = " 1, 2", x = p, y = 4*p*q^3, type = " Non-IBD "),
					data.table(gt = " 2, 2", x = p, y = q^4, type = " Non-IBD "),
					data.table(gt = " 0, 0", x = p, y = p^3, type = "IBD"),
					data.table(gt = " 0, 1", x = p, y = 2*q*p^2, type = "IBD"),
					data.table(gt = " 0, 2", x = p, y = 0, type = "IBD"),
					data.table(gt = " 1, 1", x = p, y = p^2*q + p*q^2, type = "IBD"),
					data.table(gt = " 1, 2", x = p, y = 2*p*q^2, type = "IBD"),
					data.table(gt = " 2, 2", x = p, y = q^3, type = "IBD"))

gg = ggplot(e) +
	facet_grid(.~type) +
	geom_line(aes(x, y, colour = gt), size = 1) +
	scale_x_continuous(expand = c(0.025, 0)) +
	scale_y_continuous(expand = c(0.0125, 0)) +
	scale_color_wsj() +
	theme_few() +
	theme(aspect.ratio = 2,
				strip.text = element_text(face = "bold"),
				panel.grid = element_line(colour = "grey", size = 0.25),
				panel.grid.major = element_line(colour = "grey", size = 0.25),
				panel.grid.minor = element_blank(),
				panel.margin.x = unit(1, "cm"),
				legend.title = element_blank()) +
	xlab("Allele frequency") + ylab("Probability")

ggsave(gg, filename = "__plot.expected.pairprop.pdf", height = 8, width = 10)





