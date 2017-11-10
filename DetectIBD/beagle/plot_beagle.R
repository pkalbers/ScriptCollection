

library(data.table)
library(ggplot2)
library(ggthemes)
library(colorspace)

se <- function(x) sqrt(var(x)/length(x))


load("match.beagle.all_sub_segments.RData")

H = subH
P = subP


d = rbind(data.table(side="LHS", fk = H$fk, type = "(a) True haplotypes ",   det = H$position - H$b.lhs + 1,   tru = H$position - H$t.lhs + 1),
					data.table(side="LHS", fk = P$fk, type = "(b) Phased haplotypes ", det = P$position - P$b.lhs + 1,   tru = P$position - P$t.lhs + 1),
					data.table(side="RHS", fk = H$fk, type = "(a) True haplotypes ",   det = H$b.rhs - H$position + 1,   tru = H$t.rhs - H$position + 1),
					data.table(side="RHS", fk = P$fk, type = "(b) Phased haplotypes ", det = P$b.rhs - P$position + 1,   tru = P$t.rhs - P$position + 1))

# d = split(d, list(d$type))
# k = sapply(d, function(x) unique(x$idx))
# k = Reduce(intersect, k)




### STATS
se <- function(x) sqrt(var(x)/length(x))

rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}


# over/underest
Z = P
x = 1:nrow(Z)
x = which(Z$fk == 25)
ou = data.table(r = round(Z$b.rhs[x] / 10), tr = round(Z$t.rhs[x] / 10),
								l = round(Z$b.lhs[x] / 10), tl = round(Z$t.lhs[x] / 10))

r = which(ou$r > ou$tr)
l = which(ou$l < ou$tl)
(length(l)+length(r)) / (nrow(Z[x])*2) * 100  #  H : 47.03292 (45.67136, 46.5345)  P : 45.9847 (37.25316, 46.21281)

r = which(ou$r < ou$tr)
l = which(ou$l > ou$tl)
(length(l)+length(r)) / (nrow(Z[x])*2) * 100  #  H : 49.8754 (48.35232, 51.07528)  P : 51.06569 (58.75466, 51.42911)

r = which(ou$r == ou$tr)
l = which(ou$l == ou$tl)
(length(l)+length(r)) / (nrow(Z[x])*2) * 100  #  H : 3.091682  P : 2.949609



x = by(abs(d$tru - d$det)/1e6, list(d$type, d$fk), mean)
array(round(x, 3), dim(x), dimnames(x))
# 2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25
# (a) True haplotypes    1.381 0.786 0.549 0.446 0.382 0.337 0.296 0.271 0.255 0.231 0.225 0.226 0.199 0.186 0.199 0.190 0.196 0.175 0.189 0.171 0.162 0.166 0.152 0.161
# (b) Phased haplotypes  2.301 1.046 0.689 0.533 0.433 0.383 0.339 0.302 0.274 0.245 0.238 0.239 0.213 0.192 0.213 0.196 0.212 0.193 0.203 0.175 0.169 0.182 0.159 0.171

x = by(abs(d$tru - d$det)/1e6, list(d$type, d$fk), se)
array(round(x, 3), dim(x), dimnames(x))
# 2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25
# (a) True haplotypes    0.020 0.008 0.005 0.003 0.003 0.002 0.002 0.002 0.002 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
# (b) Phased haplotypes  0.063 0.023 0.013 0.010 0.006 0.006 0.005 0.004 0.004 0.003 0.003 0.003 0.003 0.002 0.003 0.003 0.003 0.002 0.002 0.002 0.002 0.002 0.002 0.002

by(d, d$type, function(x) cor(x$tru, x$det, method = "p")^2)
by(d, list(d$type, d$side), function(x) cor(x$tru, x$det, method = "p")^2)
by(d, list(d$type, d$side), function(x) cor(x$tru, x$det, method = "s"))

by(d, d$type, function(x) { q = table(x$tru < x$det); (q/sum(q))*100 })
by(d, c(d$side, d$type), function(x) { q = table(x$tru < x$det); (q/sum(q))*100 })


x = by(d, list(d$type), function(x) cor(x$tru, x$det, method = "p")^2)
t(array(x, dim(x), dimnames(x)))

x = by(d, list(d$type, d$fk), function(x) { cor(x$tru, x$det, method = "p")^2 })
x = t(array(x, dim(x), dimnames(x)))

x = by(q, list(q$type, q$fk), function(x) { cor(x$tru, x$det, method = "p")^2 })
x = t(array(x, dim(x), dimnames(x)))


y = by(d, list(d$type), function(x) rmsle(x$tru, x$det))
t(array(y, dim(y), dimnames(y)))

y = by(d, list(d$type, d$fk), function(x) { rmsle(x$tru, x$det) })
y = t(array(y, dim(y), dimnames(y)))

y = by(q, list(q$type, q$fk), function(x) { rmsle(x$tru, x$det) })
y = t(array(y, dim(y), dimnames(y)))

cbind(x, y)

z = by(d, list(d$type, d$fk), function(x) nrow(x)/nr*100)
z = t(array(z, dim(z), dimnames(z)))




### hist, one-sided

d$map = (d$det / d$tru)

p = d

# del = which(p$map > 10)
# if (length(del) > 0) p = p[-del,]

k = which(p$fk %in% c(2, 5, 10, 15, 20, 25))
p = p[k, ]

p = split(p, list(p$type, p$fk))
p = lapply(p, function(z) {
	x = seq(0, 10, by=0.05)
	y = ecdf(z$map)(x)
	data.table(type = z$type[1], fk = z$fk[1], x=x, y=y)
})
p = rbindlist(p)

p$fk = factor(p$fk, levels = c(2:25), ordered = T)

hist = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	geom_vline(xintercept=c(1), colour="white", size = 1.25) +
	geom_line(aes(x, y, color = fk), size = 1, alpha = 0.9) +
	geom_vline(xintercept=c(1), colour="black", size = 0.75, linetype = "22") +
	scale_y_continuous(breaks = seq(0, 1, by=0.2), minor_breaks = seq(0, 1, by=0.1)) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
	coord_cartesian(xlim = c(-0.005, 2.7505), ylim = c(-0.005, 1.005), expand = F) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey90", "purple"))(3)[2],
																										 colorRampPalette(c("grey70", "navy"))(3)[2],
																										 colorRampPalette(c("grey20", "red"))(3)[2],
																										 colorRampPalette(c("grey95", "goldenrod4"))(3)[2])))(6)) +
	theme_few() +
	theme(aspect.ratio=1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
				panel.grid.minor.y = element_line(size = 0.5, colour = "grey80"),
				panel.grid.major.x = element_line(size = 0.5, colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance to true breakpoint") +
	ylab("CDF") +
	guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
hist

ggsave(filename = "____plot.beagle.hist.pdf", plot = hist, width = 9, height = 6.1)




### scatter

brk = seq(log10(1), log10(100e6), length.out = 151)

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


tm = rbind(data.table(x = c(tik*-1 -1,tik+1), y = 0, xend = c(tik*-1 -1,tik+1), yend = 5),
					 data.table(y = tik, x = 150, yend = c(tik,tik), xend = 145),
					 data.table(y = tik, x = -150, yend = c(tik,tik), xend = -145))


scatter = ggplot(data = p) +
	facet_wrap(~type, ncol = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	geom_abline(slope = c(-1, 1), intercept = c(0,0), colour="black", alpha = 0.25) +
	#scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 2.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	#scale_fill_gradient2(low = "grey", mid = "darkblue", high = "orange", na.value = "grey85", midpoint = 1.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_fill_gradientn(colours = c("white", "snow", "yellow", "gold", "orange", "orangered1", "orangered3", "firebrick3", "firebrick4"), na.value = "grey85", limits = c(0, 5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 200/401,
				panel.background = element_rect(fill = "grey60"),
				legend.background = element_rect(fill = "grey85"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				legend.title = element_blank()) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")
scatter

ggsave(filename = "____plot.beagle.scat.pdf", plot = scatter, width = 9, height = 6.1)






###
### length
###

library(data.table)

load("../data.OutOfAfricaHapMap20.RData")

#M = fread("../ooa.marker.txt", header = T, stringsAsFactors = T)
M = fread("../ooa.marker.txt", header = T, stringsAsFactors = T)

GEN = M$GenDist


load("match.beagle.all_sub_segments.RData")

H = subH
P = subP


del = which(H$b.lhs < POS[2] | H$b.rhs > tail(POS, n = 2)[1] | H$wall)  ##  35474
if (length(del) > 0) H = H[-del, ]

del = which(P$b.lhs < POS[2] | P$b.rhs > tail(POS, n = 2)[1] | P$wall)  ##  8892
if (length(del) > 0) P = P[-del, ]



pos = round(POS)

i = anyDuplicated(pos)
while (i != 0) {
	if (i == 1) {
		pos[i] = pos[i] - 1
	} else {
		pos[i] = pos[i] + 1
	}
	i = anyDuplicated(pos)
}



h.phy = ((H$b.rhs - H$b.lhs)+1) / 1e6
p.phy = ((P$b.rhs - P$b.lhs)+1) / 1e6
t.phy = ((H$t.rhs - H$t.lhs)+1) / 1e6
f.phy = ((P$t.rhs - P$t.lhs)+1) / 1e6

h.gen = (GEN[ match((H$b.rhs), pos) ] - GEN[ match((H$b.lhs), pos) ]) + 1e-08
p.gen = (GEN[ match((P$b.rhs), pos) ] - GEN[ match((P$b.lhs), pos) ]) + 1e-08
t.gen = (GEN[ match((H$t.rhs), pos) ] - GEN[ match((H$t.lhs), pos) ]) + 1e-08
f.gen = (GEN[ match((P$t.rhs), pos) ] - GEN[ match((P$t.lhs), pos) ]) + 1e-08


q = rbind(data.table(fk = H$fk, type = "(a) True haplotypes",    x = h.phy, mode = "A. Physical length"),
					data.table(fk = H$fk, type = "(a) True haplotypes",    x = h.gen, mode = "B. Genetic length"),
					data.table(fk = P$fk, type = "(b) Phased haplotypes",  x = p.phy, mode = "A. Physical length"),
					data.table(fk = P$fk, type = "(b) Phased haplotypes",  x = p.gen, mode = "B. Genetic length"),
					data.table(fk = H$fk, type = "True IBD",                x = t.phy, mode = "A. Physical length"),
					data.table(fk = H$fk, type = "True IBD",                x = t.gen, mode = "B. Genetic length"))#,
					# data.table(fk = P$fk, type = "True IBD x",                x = f.phy, mode = "A. Physical length"),
					# data.table(fk = P$fk, type = "True IBD x",                x = f.gen, mode = "B. Genetic length"))


# q = split(q, list(q$type, q$mode))
# cor(q$`True IBD.A. Physical length`$x, q$`(a) FGT, true haplotypes.A. Physical length`$x, method = 'p')^2
# cor(q$`True IBD.A. Physical length`$x, q$`(a) FGT, true haplotypes.A. Physical length`$x, method = 'p')^2


x = by(q$x, list(q$type, q$mode), median)
array(x, dim(x), dimnames(x))

x = by(q$x, list(q$fk, q$type, q$mode), median)
array(x, dim(x), dimnames(x))



q = split(q, list(q$fk, q$type, q$mode))

# sub = sample(match$index, 10000)
# q = lapply(q, function(x) {
# 	k = which(x$idx %in% sub)
# 	x[k, ]
# })
# p = rbindlist(q)

q = lapply(q, function(x) {
	z = x
	w = fivenum(x$x)
	z$low = w[2]
	z$med = w[3]
	z$upp = w[4]
	z[1,]
})
q = rbindlist(q)


#p$fk = factor(p$fk, levels = (2:25), ordered = T)
#q$fk = factor(q$fk, levels = (2:25), ordered = T)

tl = c("True IBD", "(a) True haplotypes","(b) Phased haplotypes")

q$type = factor(q$type, levels = tl, ordered = T)

ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))

len = ggplot(q) + #[sample(1:nrow(p), 10000),]) +
	facet_wrap(~mode, scales = "free", ncol = 1) +
	#geom_violin(aes(fk, x, fill = type), size = 0, alpha = 0.9) +
	#geom_pointrange(data=q, aes(x=fk, y=med, ymin=low, ymax=upp, group=type), color = "black", size = 2/3, shape='|', position=position_dodge(width = 0.9)) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=type), size = 1.5, position=position_dodge(width = 0.8)) +
	geom_point(aes(x=fk, y=med, color=type), size = 3, shape=18, position=position_dodge(width = 0.8)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), minor_breaks = mb) +
	scale_x_continuous(breaks = 2:25) +
	scale_color_manual(values = c("grey20", "royalblue1", "purple", "limegreen")) +
	#scale_fill_manual(values = c("turquoise3", "purple", "limegreen", "grey50")) +
	coord_cartesian(ylim = c(0.025, 12.5), xlim = c(1.5, 25.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "vertical",
				legend.text = element_text(size = 9),
				legend.key.height = unit(0.1, "cm"),
				legend.background = element_rect(fill = "white", colour = "grey50", size = 0.5),
				#panel.spacing.x = unit(-1, "points"),
				strip.text = element_text(face = "bold", hjust = 0),
				#strip.text.y = element_blank(),
				#axis.ticks.y = element_blank(),
				panel.grid = element_line(colour = "grey80", size = 1/3),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") +
	ylab("Genetic length (cM)                                     Physical length (Mb)") +
	guides(color = guide_legend(override.aes = list(size = 5)))
len

ggsave(len, filename = "____plot.beagle.length.pdf", height=7, width = 9)



##   stats



x = by(q$x, list(q$type, q$mode), median)
array(x, dim(x), dimnames(x))




p = q
p = split(p, list(p$type, p$fk))

for (tag in names(p)) cat(tag, lm(log(x) ~ fk, p[[tag]])$coefficients, "\n")


