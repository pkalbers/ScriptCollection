

expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}

expected.age = function(k, n) {
	if (length(k) > 1) return(sapply(k, expected.age, n))
	j = 2:n
	s = sum( choose(n - j, k - 1) * ((n - j + 1) / (n * (j - 1))) )
	2 * choose(n - 1, k)^-1 * s
}


get.trans.prob = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	t = expected.age(k, n)
	r = d * rec.rate 
	p = (4/6) * (1 - exp(-6 * t))
	
	q11 = exp(-4 * N * r * t)
	q10 = 1 - q11
	q01 = (p / (1 - p)) * q10
	q00 = 1 - q01
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("ibd", "non"), c("ibd", "non")))
}


dist = c(0, 10^seq(log10(1), log10(1e5), length.out = 101))
freq = c(2, 5, 10, 20, 30, 40, 50, 250)


# plot trans prob
g = NULL
n = 5000
for (d in dist) {
	for (k in freq) {
		tmp = as.data.table(melt(get.trans.prob(d, k, n)))
		names(tmp) = c("from", "to", "prob")
		tmp$d = d
		tmp$k = k
		g = rbind(g, tmp)
	}
}

g$tag = sprintf("bolditalic(%s) %%->%% bolditalic(%s)", g$from, g$to)
g$frq = sprintf("%s %%", as.character(g$k / n * 100))

g = g[-(which(g$from == "non")), ]


gg = ggplot(g) + 
	facet_wrap(~tag, ncol = 2, scales = "free_x", labeller = label_parsed) + 
	geom_line(aes(x=d, y=prob, colour = factor(frq)), size = 2/3) +
	scale_x_continuous(expand = c(0.02, 0), breaks = seq(0, 1e5, by = 10000), labels = seq(0, 1e5, by = 10000) * 1e-8 * 100) +
	scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 1, by = 0.1)) +
	scale_color_brewer(palette = "Dark2") +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				#axis.text = element_text(size = 7) ,
				axis.title.x = element_text(margin=margin(10,0,0,0)),
				axis.title.y = element_text(margin=margin(0,10,0,0)),
				panel.margin = unit(0.5, "cm"),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank(),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold")) +
	xlab("Genetic distance (cM)") + 
	ylab("Transition probability")

gg

ggsave(gg, filename = "_plot.transmatrix.pdf", width = 9, height = 4.5)


