

library(ggplot2)
library(ggthemes)
library(data.table)





eij = function(f, g, i, j, Y, N) {
	p = NULL
	
	if (j == 0) {
		if (i == 0) p = (1-Y) * f  #  (((1-em)*((1-e0)^2)) + em*(e0^2)) * f
		if (i == 1) p = Y * f
		if (i == 2) p = 0 * f  #  ((1-em)*(e0^2) + (em * ((1-e0)^2))) * f
	}
	if (j == 1) {
		if (i == 0) p = (0.5 * N) * f  #  (((1-em)*(e1*(1-e0))) + (em * (e0*(1-e1)))) * f
		if (i == 1) p = (1 - N) * f
		if (i == 2) p = (0.5 * N) * f  #  (((1-em) * (e0*(1-e1))) + (em * (e1*(1-e0)))) * f
	}
	if (j == 2) {
		if (i == 0) p = 0 * f  #  ((1-em)*(e1^2) + (em * ((1-e1)^2))) * f
		if (i == 1) p = Y * f
		if (i == 2) p = (1-Y) * f  #  ((1-em)*((1-e1)^2) + (em * (e1^2))) * f
	}
	
	p
}

eij = function(f, g, i, j, e0, e1, em=0.01) {
	p = NULL
	
	if (j == 0) {
		if (i == 0) p = (1-e0)^2 * f  #  (((1-em)*((1-e0)^2)) + em*(e0^2)) * f
		if (i == 1) p = 2*e0*(1-e0) * f
		if (i == 2) p = (e0^2) * f  #  ((1-em)*(e0^2) + (em * ((1-e0)^2))) * f
	}
	if (j == 1) {
		if (i == 0) p = (e1*(1-e0)) * f  #  (((1-em)*(e1*(1-e0))) + (em * (e0*(1-e1)))) * f
		if (i == 1) p = (e0*e1 + (1-e0)*(1-e1)) * f
		if (i == 2) p = (e0*(1-e1)) * f  #  (((1-em) * (e0*(1-e1))) + (em * (e1*(1-e0)))) * f
	}
	if (j == 2) {
		if (i == 0) p = (e1^2) * f  #  ((1-em)*(e1^2) + (em * ((1-e1)^2))) * f
		if (i == 1) p = 2*e1*(1-e1) * f
		if (i == 2) p = ((1-e1)^2) * f  #  ((1-em)*((1-e1)^2) + (em * (e1^2))) * f
	}
	
	p
}

eij = function(f, g, i, j, e0, e1, em=0.1) {
	p = NULL
	
	# if (j == 0) {
	# 	if (i == 0) p = ( (1-em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) +   ((em)^2)*(e0^2) ) * f
	# 	if (i == 1) p =   ( (em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) +   ((em)^2)*(e0^2) ) * f
	# 	if (i == 2) p =   ( (em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) + ((1-em)^2)*(e0^2) ) * f
	# }
	if (j == 0) {
		if (i == 0) p = ( (1-em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) +   ((em)^2)*(e0^2) ) * f
		if (i == 1) p =   ( (em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) +   ((em)^2)*(e0^2) ) * f
		if (i == 2) p =   ( (em)^2*((1-e0)^2) + (2*(1-em)*em)*2*e0*(1-e0) + ((1-em)^2)*(e0^2) ) * f
	}
	if (j == 1) {
		if (i == 0) p = (e1*(1-e0)) * f  #  (((1-em)*(e1*(1-e0))) + (em * (e0*(1-e1)))) * f
		if (i == 1) p = (e0*e1 + (1-e0)*(1-e1)) * f
		if (i == 2) p = (e0*(1-e1)) * f  #  (((1-em) * (e0*(1-e1))) + (em * (e1*(1-e0)))) * f
	}
	if (j == 2) {
		if (i == 0) p = (e1^2) * f  #  ((1-em)*(e1^2) + (em * ((1-e1)^2))) * f
		if (i == 1) p = 2*e1*(1-e1) * f
		if (i == 2) p = ((1-e1)^2) * f  #  ((1-em)*((1-e1)^2) + (em * (e1^2))) * f
	}
	
	p
}



pd = function(e0, e1, em = 0.05) {
	q = seq(0,1, by=0.001)
	p = 1 - q
	
	g0 = p^2
	g1 = 2*p*q
	g2 = q^2
	
	
	f00 = eij(g0, g0, 0, 0, e0, e1, em)
	f10 = eij(g1, g0, 1, 0, e0, e1, em)
	f20 = eij(g2, g0, 2, 0, e0, e1, em)
	sum = f00 + f10 + f20
	
	f00 = f00 / sum
	f10 = f10 / sum
	f20 = f20 / sum
	
	
	f01 = eij(g0, g1, 0, 1, e0, e1, em)
	f11 = eij(g1, g1, 1, 1, e0, e1, em)
	f21 = eij(g2, g1, 2, 1, e0, e1, em)
	sum = f01 + f11 + f21
	
	f01 = f01 / sum
	f11 = f11 / sum
	f21 = f21 / sum
	
	
	f02 = eij(g0, g2, 0, 2, e0, e1, em)
	f12 = eij(g1, g2, 1, 2, e0, e1, em)
	f22 = eij(g2, g2, 2, 2, e0, e1, em)
	sum = f02 + f12 + f22
	
	f02 = f02 / sum
	f12 = f12 / sum
	f22 = f22 / sum
	
	tg = function(gt) sprintf("bold(g[%s])", gt)
	og = function(gt) sprintf("bold(tilde(g)[%s])", gt)
	
	d = rbind(data.table(frq = q, obs = og("0"), tru = tg("0"), p = f00),
						data.table(frq = q, obs = og("1"), tru = tg("0"), p = f10),
						data.table(frq = q, obs = og("2"), tru = tg("0"), p = f20),
						data.table(frq = q, obs = og("0"), tru = tg("1"), p = f01),
						data.table(frq = q, obs = og("1"), tru = tg("1"), p = f11),
						data.table(frq = q, obs = og("2"), tru = tg("1"), p = f21),
						data.table(frq = q, obs = og("0"), tru = tg("2"), p = f02),
						data.table(frq = q, obs = og("1"), tru = tg("2"), p = f12),
						data.table(frq = q, obs = og("2"), tru = tg("2"), p = f22))
	
	d
}


a = pd(1/1000, 1/1000)
b = pd(5/1000, 10/1000)
c = pd(10/1000, 5/1000)

sa = sprintf("     = %.3f,      = %.3f  ", 1/1000, 1/1000)
sb = sprintf("     = %.3f,      = %.3f  ", 5/1000, 10/1000)
sc = sprintf("     = %.3f,      = %.3f  ", 10/1000, 5/1000)

d = rbind(cbind(type="a", a), cbind(type="b", b), cbind(type="c", c))

d$type = factor(d$type)


gg = ggplot(d) + 
	facet_grid(obs~tru, labeller = label_parsed) + 
	geom_line(aes(x = frq, y = p, linetype = type, colour = type), size = 0.6) +
	scale_linetype_manual(values = c("solid", "31", "11"), labels = c(a = sa, b = sb, c = sc)) +
	scale_colour_manual(values = c("red", "forestgreen", "blue"), labels = c(a = sa, b = sb, c = sc)) +
	scale_x_continuous(expand = c(0.04, 0)/2, breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_y_continuous(expand = c(0.04, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "", "0.5", "", "1")) +
	theme_few()+
	theme(#aspect.ratio = 1/3,
		legend.title = element_blank(),
		legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.width = unit(1.5, "cm"),
		legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.margin = unit(0.1, "cm"),
		strip.text.x = element_text(face = "bold", angle = 0),
		strip.text.y = element_text(face = "bold", angle = 0),
		panel.grid = element_line(colour = "grey85", linetype = "solid", size = 0.25),
		panel.grid.major = element_line(colour = "grey85", linetype = "solid", size = 0.25),
		panel.grid.minor = element_blank(),
		panel.margin.y = unit(0.5, "mm"),
		panel.border = element_rect(fill = NA, color = "black", size = 0.5),
		axis.text = element_text(size = 8),
		plot.margin = unit(c(1, 1, 1, 1)/2, "mm")) +
	xlab("Allele frequency") +
	ylab("Expected proportion")
gg

ggsave(gg, filename = "_plot.errmod_genotype.pdf", width = 9, height = 4.5)
ggsave(gg, filename = "_plot.errmod_allele.pdf", width = 9, height = 4.5)



a = pd(1/1000, 1/1000, 0.01)
b = pd(1/1000, 1/1000, 0.05)
c = pd(1/1000, 1/1000, 0.10)

sa = sprintf("     = %.2f  ", 0.01)
sb = sprintf("     = %.2f  ", 0.05)
sc = sprintf("     = %.2f  ", 0.10)

d = rbind(cbind(type="a", a), cbind(type="b", b), cbind(type="c", c))

d$type = factor(d$type)

gg = ggplot(d) + 
	facet_grid(obs~tru, labeller = label_parsed) + 
	geom_line(aes(x = frq, y = p, linetype = type), size = 0.6) +
	scale_linetype_manual(values = c("solid", "31", "11"), labels = c(a = sa, b = sb, c = sc)) +
	scale_x_continuous(expand = c(0.04, 0)/2, breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_y_continuous(expand = c(0.04, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "", "0.5", "", "1")) +
	theme_few()+
	theme(#aspect.ratio = 1/3,
		legend.title = element_blank(),
		legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.width = unit(1.5, "cm"),
		legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.margin = unit(0.1, "cm"),
		strip.text.x = element_text(face = "bold", angle = 0),
		strip.text.y = element_text(face = "bold", angle = 0),
		panel.grid = element_line(colour = "grey85", linetype = "solid", size = 0.25),
		panel.grid.major = element_line(colour = "grey85", linetype = "solid", size = 0.25),
		panel.grid.minor = element_blank(),
		panel.margin.y = unit(0.5, "mm"),
		panel.border = element_rect(fill = NA, color = "black", size = 0.5),
		axis.text = element_text(size = 8),
		plot.margin = unit(c(1, 1, 1, 1)/2, "mm")) +
	xlab("Allele frequency") +
	ylab("Expected proportion")
gg

ggsave(gg, filename = "_plot.errmod_allele_ext.pdf", width = 9, height = 4.5)



z = split(d, d$type)
z = lapply(z, function(d) {
	x = split(d, d$tru)
	z = list()
	for (tag in names(x)) {
		f = c()
		if (grepl("0", tag)) f = (1 - x[[tag]]$frq)^2
		if (grepl("1", tag)) f = 2 * (x[[tag]]$frq) * (1 - x[[tag]]$frq)
		if (grepl("2", tag)) f = (x[[tag]]$frq)^2
		y = lapply(x[tag], function(d, f) {
			y = sapply(split(x[[tag]], x[[tag]]$obs), function(d) sum(d$p * f))
			y / sum(y)
		}, f)
		z[[tag]] = sprintf("%.4f", unlist(y, use.names = F) * 100)
	}
	z
})
z


e0 = 10/1000
em = 0.01

(1-e0)^2
2*e0*(1-e0)
(e0^2)

((((1-e0)^2)) + em*((e0)^2))
2*e0*(1-e0)
((e0^2) + (em * ((1-e0)^2)))


