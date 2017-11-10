

library(data.table)
library(ggplot2)
library(ggthemes)


Ne = 2 * 10000

marker = fread("./truH.marker.txt", header = T)
times = fread("./OutOfAfricaHapMap20.times.txt", header = T)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}

se <- function(x) sqrt(var(x)/length(x))



load("results.data.RData")

d = cbind(d, times[d$MarkerID + 1, ])

d$clock = ""
d$method = ""

d$clock[which(grepl("^.+_1_0$", d$tag))] = "Mutation clock"
d$clock[which(grepl("^.+_0_1$", d$tag))] = "Recombination clock"
d$clock[which(grepl("^.+_1_1$", d$tag))] = "Combined clock"

d$method[which(grepl("^fgt_(err|tru)H_.+$", d$tag))] = "FGT"
d$method[which(grepl("^fgt_(err|tru)P_.+$", d$tag))] = "FGT*"
d$method[which(grepl("^dgt_.+$", d$tag))] = "DGT"
d$method[which(grepl("^hmm_.+$", d$tag))] = "HMM"



p = d
p = p[which(p$error == T),]

load("result.truth.RData")

p = rbind(p, sim)


p = p[order(p$clock, p$method),]
p = split(p, list(p$clock, p$method))
p = rbindlist(p)
p = as.data.frame(as.list(p))

p$clock = factor(p$clock, levels = c("Mutation clock", "Recombination clock", "Combined clock" ), ordered = T)
p$method = factor(p$method, levels = c("True IBD", "FGT", "FGT*", "DGT", "HMM"), ordered = T)




del = which(p$Lower > p$Upper)  ; if (length(del) > 0) p=p[-del,]
del = which(p$Robust > p$Upper)  ; if (length(del) > 0) p=p[-del,]
del = which(p$Robust < p$Lower)  ; if (length(del) > 0) p=p[-del,]
del = which(p$PostMode <= 1e-07)  ; if (length(del) > 0) p=p[-del,]
table(p$method)

k = lapply(split(p$MarkerID, list(p$method, p$clock)), unique)
k = Reduce(intersect, k)
tmp = k
key = intersect(tmp, k)
save(key, file = "tmp.RData")

p$mid = exp(log(p$node.time) + ((log(p$parent.time) - log(p$node.time)) / 2))
p$x = p$mid # p$node.time

pp = p[which(p$node.time >= 1 & p$node.time <= 1000), ]



mx = c(0.5, 20000)
tr = rbind(data.table(method = "True IBD", clock = "Mutation clock",      x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]),
					 data.table(method = "True IBD", clock = "Recombination clock", x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]),
					 data.table(method = "True IBD", clock = "Combined clock",      x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]))
tr$clock = factor(tr$clock, levels = c("Mutation clock", "Recombination clock", "Combined clock" ), ordered = T)
tr$method = factor(tr$method, levels = c("True IBD", "FGT", "FGT*", "DGT", "HMM"), ordered = T)


tm = rbind(data.table(x = 10^(0:6), y = 0.0001, xend = 10^(0:6), yend = 0.7),
					 data.table(y = 10^(0:6), x = 0.0001, yend = 10^(0:6), xend = 0.7))


scat = ggplot(p) +
	facet_grid(method~clock) +
	geom_raster(aes(x, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(10/100, 10/100)) +
	geom_smooth(data = pp, aes(x, node.time), method = "lm", size = 1/2, color = "grey30", fullrange = T) +
	#geom_smooth(data = pp, aes(x, mid), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	geom_smooth(data = pp, aes(x, parent.time), method = "lm", size = 1/2, color = "grey30", fullrange = T) +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = T) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", fullrange = T) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy", "black"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "grey50", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,1), legend.position=c(1-0.99,0.995),
				legend.background = element_rect(fill = NA, colour = NA, size = NA),
				#legend.position = "bottom",
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age")
scat


ggsave(scat, filename = "_plot.scat.tru.pdf", height = 11, width = 7.2)
ggsave(scat, filename = "_plot.scat.err.pdf", height = 11, width = 7.2)












a = sapply(split(p, list(p$clock, p$method)), function(x) {
	data.table("Child node"  = cor(x$PostMode * Ne, x$node.time, method = "s"),
						 "Midpoint"    = cor(x$PostMode * Ne, x$mid, method = "s"),
						 "Parent node" = cor(x$PostMode * Ne, x$parent.time, method = "s"))
	
})
a = t(array(a, dim = dim(a), dimnames = dimnames(a)))
a

b = sapply(split(p, list(p$clock, p$method)), function(x) {
	data.table("Child node"  = rmsle(x$PostMode * Ne, x$node.time),
						 "Midpoint"    = rmsle(x$PostMode * Ne, x$mid),
						 "Parent node" = rmsle(x$PostMode * Ne, x$parent.time))
	
})
b = t(array(b, dim = dim(b), dimnames = dimnames(b)))
b


x = by(p, list(p$Fk, p$clock, p$method), function(x) cor(x$PostMode * Ne, x$node.time, method = "s"))
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(p, list(p$Fk, p$clock, p$method), function(x) rmsle(x$PostMode * Ne, x$node.time))
array(x, dim = dim(x), dimnames = dimnames(x))


x = by(p, list(p$Fk, p$clock, p$method), nrow)
array(x, dim = dim(x), dimnames = dimnames(x))



q = data.table(fk = p$Fk, clock = p$clock, method = p$method, est = p$PostMode * Ne, tru = p$mid, t0 = p$node.time, t1 = p$parent.time)

q$map = log(q$est / q$t0) / log(q$t1 / q$t0)

tmp = q
tmp = tmp[which(tmp$fk %in% c(0, 2, 5, 10, 15, 20, 25)), ]
tmp = split(tmp, list(tmp$clock, tmp$method, tmp$fk))
tmp = lapply(tmp, function(z) {
	if (nrow(z) == 0) return(NULL)
	x = seq(-10, 10, by=0.05)
	y = ecdf(z$map)(x)
	data.table(fk = z$fk[1], clock = z$clock[1], method = z$method[1], x=x, y=y)
})
q = rbindlist(tmp)

x = which(q$method == "True IBD")
r = q[x,]
r$method = NULL
q = as.data.frame(q)[-x,]


hist = ggplot(q) + 
	facet_grid(method~clock) +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	# geom_vline(xintercept = 0, size = 3/4, linetype = "21", color = "green4") +
	# geom_vline(xintercept = 1, size = 3/4, linetype = "21", color = "sienna") +
	geom_vline(xintercept = c(0,1), size = 2/4, linetype = "21") +
	#geom_vline(xintercept = 0.5, size = 3/4, color = "grey60") +
	#geom_histogram(aes(est), binwidth = 0.1, alpha = 1/2, fill = "grey20") +
	geom_line(aes(x, y, color = factor(fk)), size = 2/3) +
	geom_line(data = r, aes(x, y), color = "black", size = 1, alpha = 1/3) +
	scale_x_continuous(breaks = seq(-10, 10, by = 0.5), labels = as.character(seq(-10, 10, by = 0.5))) +
	scale_y_continuous(breaks = seq(0,1,by = 1/4)) +
	#scale_y_continuous(breaks = seq(0, 10000, by = 250), labels = format(seq(0, 10000, by = 250), big.mark = ",")) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey90", "maroon1"))(3)[2],
																										 colorRampPalette(c("grey70", "navy"))(3)[2],
																										 colorRampPalette(c("grey20", "red"))(3)[2],
																										 colorRampPalette(c("grey95", "darkorange3"))(3)[2])))(6)) +
	coord_cartesian(expand = F, xlim = c(-1.75, 2.75), ylim = c(-0.025, 1.025)) +
	theme_few() +
	theme(aspect.ratio = 1/2,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title = element_blank(),
				legend.text = element_text(size = 8),
				legend.key.height = unit(4, "mm"),
				legend.box.background = element_rect(fill = "white", size = 1/3, colour = "grey20"),
				legend.position = c(0.005, 0.99), legend.justification = c(c(0,1)), legend.margin = margin(0,1,1,1, "mm"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.y = element_blank()) + 
	xlab("Relative age estimate (log-scale)") + ylab("CDF") +
	guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
hist

ggsave(hist, filename = "_plot.hist.pdf", height = 5.5, width = 10)













q = data.table(clock = d$clock,
							 method = d$method,
							 con = d$Lower * Ne,
							 dis = d$Upper * Ne,
							 est = d$PostMode * Ne,
							 t0 = d$node.time,
							 t1 = d$parent.time,
							 err = d$error)

q$est = log(q$est / q$t0) / log(q$t1 / q$t0)


ggplot(q) +
	facet_grid(method~clock) +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	# geom_vline(xintercept = 0, size = 3/4, linetype = "21", color = "green4") +
	# geom_vline(xintercept = 1, size = 3/4, linetype = "21", color = "sienna") +
	geom_vline(xintercept = c(0,1), size = 3/4, linetype = "21") +
	geom_vline(xintercept = 0.5, size = 3/4, color = "grey60") +
	geom_histogram(data = q[which(q$err == F), ], aes(est), binwidth = 0.1, alpha = 1/2, fill = "seagreen") +
	geom_histogram(data = q[which(q$err == T), ], aes(est), binwidth = 0.1, alpha = 1/2, fill = "red3") +
	scale_x_continuous(breaks = -10:10) +
	scale_y_continuous(breaks = seq(0, 10000, by = 250), labels = format(seq(0, 10000, by = 250), big.mark = ",")) +
	coord_cartesian(expand = F, xlim = c(-3.5, 4.5), ylim = c(0, 1550)) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 9),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.x = element_line(size = 1/2, colour = "grey90")) +
	xlab("Relative age estimate (log-scale)") + ylab("Count")

ggsave(filename = "_plot.hist.pdf", height = 10, width = 8, dpi = 300)


z = by(p, list(p$clock, p$method), function(p) length(which(p$PostMode*Ne < p$node.time)) / nrow(p) * 100)
array(z, dim = dim(z), dimnames = dimnames(z))
z = by(p, list(p$clock, p$method), function(p) length(which(p$PostMode*Ne > p$parent.time)) / nrow(p) * 100)
array(z, dim = dim(z), dimnames = dimnames(z))




q = rbindlist(lapply(split(q, q$clock), function(p) {
	tmp = NULL
	for (tag in names(p)[-1]) {
		den = density(p[[tag]], from = -5, to = 6, n = 100)
		tmp = rbind(tmp, data.table(clock = p$clock[1], type = tag, x = den$x, y = den$y))
	}
	tmp
}))


ggplot(q) +
	facet_wrap(~clock, ncol = 1) +
	geom_vline(xintercept = c(0, 1)) +
	#geom_area(aes(x=x, y = y, fill = type), alpha = 1/3, position = position_identity()) +
	#geom_line(aes(x=x, y = y, color = type), alpha = 2/3, size = 1) +
	geom_bar(aes(x=x, y=y), stat = "identity") +
	#geom_hline(yintercept = 0, color = "white", size = 1) +
	#scale_colour_manual(values = c(con = "forestgreen", dis = "saddlebrown", est = "dodgerblue")) +
	#scale_fill_manual(values = c(con = "forestgreen", dis = "saddlebrown", est = "dodgerblue1")) +
	scale_x_continuous(breaks = -10:10) +
	coord_cartesian(expand = F, xlim = c(-3.5, 4.5), ylim = c(0, max(q$y + 1))) +
	theme_few() +
	theme(aspect.ratio = 1/2)



ggplot(p) +
	#facet_grid(fk ~ .) +
	#geom_histogram(aes(est), binwidth = 0.05, alpha = 0.5, fill = "deepskyblue") +
	geom_density(aes(con, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "forestgreen", fill = "green4") +
	stat_bin(aes(dis, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "saddlebrown") +
	stat_bin(aes(est, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "deepskyblue") +
	geom_vline(xintercept = c(0, 1)) +
	scale_x_continuous(breaks = -10:10) +
	coord_cartesian(xlim = c(-3, 4)) +
	theme_few()



#ggplot(p) + geom_point(aes(time, PostMode * Ne))


ggplot(p) +
	geom_smooth(aes(time, PostMode * Ne), method = "lm") +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_linerange(aes(x = time, ymax = Upper * Ne, ymin = Lower * Ne)) +
	geom_point(aes(time, PostMode * Ne)) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F)))


#d$time = times$time[ d$MarkerID + 1 ]


d$method[which(d$method == "FGT")] = "bold('(a) FGT, true haplotypes')"
d$method[which(d$method == "DGT")] = "bold('(c) DGT')"
d$method[which(d$method == "HMM")] = "bold('(d) HMM')"

d$method[which(d$type == "truP" | d$type == "errP")] = "bold('(b) FGT, phased haplotypes')"

d$clock[which(d$clock == "Mut")]    = " bolditalic(T[M]) "
d$clock[which(d$clock == "Rec")]    = " bolditalic(T[R]) "
d$clock[which(d$clock == "MutRec")] = "bolditalic(T[MR])"

d$type[which(d$type == "truH" | d$type == "truP")] = "tru"
d$type[which(d$type == "errH" | d$type == "errP")] = "err"


d$type = sprintf("bold('%s')", d$type)


p = d[which(d$type == "tru"), ]
p = d[which(d$type == "err"), ]

# subset
q = split(p, list(p$method, p$clock))
q = lapply(q, function(p) {
	if (nrow(p) < 5000) return(p)
	p[sample(1:nrow(p), 5000), ]
})
p = rbindlist(q)


tm = rbind(data.table(x = 10^(0:6), y = 0.01, xend = 10^(0:6), yend = 0.65),
					 data.table(y = 10^(0:6), x = 0.01, yend = 10^(0:6), xend = 0.65))

gg = ggplot(p) +
	facet_grid(clock~type, labeller = label_parsed) +
	# geom_rect(data = data.table(xmin = 0.2, xmax = 50000, ymin = 0.2, ymax = 50000),
	# 					aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey80") +
	geom_raster(aes(time, PostMedian * Ne), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_bin2d(aes(time, PostMean * 20000), binwidth = c(0.075, 0.075), size = 2) +
	#stat_density2d(aes(time, PostMean * 10000, fill = ..level.., alpha = ..level..), size = 0.01, bins = 50, geom = 'polygon') +
	geom_density_2d(aes(time, PostMedian * Ne), bins = 5, size = 1/3, color = "black", alpha = 0.8) +
	geom_smooth(aes(time, PostMedian * Ne), method = "lm") +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	coord_cartesian(xlim = c(0.5, 50000), ylim = c(0.5, 50000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("white", "deepskyblue", "royalblue2", "magenta4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey75", "yellow", "orange", "orangered", "darkred"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	#scale_color_gradientn(colours = c("lightblue", "deepskyblue", "royalblue", "darkblue"), na.value = "navy", limits = c(0, 4)) +
	#scale_fill_gradient(low = "green", high = "red") +
	scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				#panel.spacing.y = unit(0.05, "cm"),
				#plot.margin = unit(c(-0.1, 0.05, 0.05, 0.05), "cm"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				#panel.background = element_rect(fill = "grey80"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 8),
				strip.text.y = element_text(face = "bold.italic", angle = 0),
				legend.justification=c(1,0), legend.position=c(0.995,1-0.995), legend.box.background = element_rect(colour = "black"),
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age") #+ labs(fill = "Density")
gg

ggsave(gg, filename = "_plot.age.check.HMM.postmean.pdf", width = 6, height = 8)
ggsave(gg, filename = "_plot.age.check.HMM.postmode.pdf", width = 6, height = 8)
ggsave(gg, filename = "_plot.age.check.HMM.postmedi.pdf", width = 6, height = 8)

ggsave(gg, filename = "_plot.age.sites.africa.tru.pdf", width = 8, height = 6.5)
ggsave(gg, filename = "_plot.age.sites.africa.err.pdf", width = 8, height = 6.5)



q = split(d, list(d$type, d$method, d$clock, d$ndiscord))
q = lapply(q, function(p) {
	if (nrow(p) < 5000) return(p)
	p[sample(1:nrow(p), 5000), ]
})
d = rbindlist(q)


x = by(d, list(d$Fk, d$type, d$method, d$clock), function(x) cor(x$time, x$PostMean, method = "s"))
array(x, dim = dim(x), dimnames = dimnames(x))

y = by(d, list(d$Fk, d$type, d$method, d$clock), function(x) rmsle(x$time, x$PostMedian * Ne))
array(y, dim = dim(y), dimnames = dimnames(y))

z = by(d, list(d$type, d$method, d$clock), function(x) mean(abs(x$time - x$PostMean * Ne)))
array(z, dim = dim(z), dimnames = dimnames(z))





x = by(d, list(d$type, d$method, d$clock), nrow)
array(x, dim = dim(x), dimnames = dimnames(x))



ggplot(p) +
	#facet_grid(clock~., labeller = label_parsed) +
	stat_summary(aes(Fk, PostMean * 20000), fun.data = "median_hilow", geom = "smooth", alpha = 0.5) +
	stat_summary(aes(Fk, time), fun.data = "median_hilow", geom = "smooth", color = "black", alpha = 0.5) +
	scale_y_log10() +
	scale_x_log10() +
	theme_few()




###
### ration of con/dis lengths
###

load("data.result.age.pairs.check.RData")

markers = fread("vanilla.marker.txt", header = T)

d$len = markers$Position[d$SegmentRHS + 1] - markers$Position[d$SegmentLHS + 1] + 1

d$Shared = as.character(d$Shared)

ggplot(d) +
	#facet_grid(clock~type) +
	geom_histogram(aes(len, fill = Shared), bins = 100, position=position_dodge()) +
	scale_x_log10() +
	scale_y_log10() +
	theme_bw() +
	xlab("Physical length") + ylab("Count")




###












markers = fread("vanilla.marker.txt", header = T, stringsAsFactors = F)

files = dir(pattern = ".+\\.age\\.sites\\..+")



file = "test.age.sites.FGT.MutRec.txt"
file = "test.age.sites.FGT.Mut.txt"
file = "test.age.sites.FGT.Rec.txt"

file = "test.age.sites.DGT.MutRec.txt"
file = "test.age.sites.DGT.Mut.txt"
file = "test.age.sites.DGT.Rec.txt"

file = "test.age.sites.HMM.MutRec.txt"
file = "test.age.sites.HMM.Mut.txt"
file = "test.age.sites.HMM.Rec.txt"

age = read.age(file, M, R)

ggplot(age) +
	geom_point(aes(time, PostMean * 10000, color = (Fk))) +
	geom_abline(intercept = c(0,0), slope = 1) +
	coord_cartesian(xlim = c(0.2, 200000), ylim = c(0.2, 200000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = format(10^(0:6), big.mark = ',', scientific = F)) +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), big.mark = ',', scientific = F))

cor(age$time, age$PostMean, method = "s")
cor(age$Fk, age$PostMean, method = "s")
cor(age$time, age$Fk, method = "s")

cor(age$time, age$PostMedian, method = "s")
cor(age$time, age$PostMode, method = "s")
cor(age$time, age$Robust, method = "s")
cor(age$time, age$Lower, method = "s")

plot(sort(age$PostMean*20000), ylim=c(0, 20000))
plot(sort(age$time), ylim=c(0,20000))

cor(age$PostMean, age$Fk)^2

plot(age$Fk, age$PostMean)


rmsle = function(a, b) {
	n = length(a)
	a = log(a + 1)
	b = log(b + 1)
	sqrt(sum((a - b)^2) / n)
}

rmsle(age$time, age$PostMean)


d = NULL

for (file in files) {
	
	method = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\1", file)
	clock  = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\2", file)
	
	phased = grepl(".+P\\.age\\..+", file)
	
	if (phased)
		method = sprintf("%s phased", method)
	
	#if (clock == "MutRec") clock = " MutRec"
	
	cat(".")
	
	age = read.age(file, M, R)
	
	age$method = method
	age$clock = clock
	
	d = rbind(d, age)
}


del = which(d$N_Others < d$N_Shared)
if (length(del) > 0)
	d = d[-del, ]


m = match(d$MarkerID, markers$MarkerID)
d$a = markers$AlleleCount1[m]
d$g = markers$GenotypeCount1[m]
del = which(d$a != d$g)
if (length(del) > 0)
	d = d[-del, ]


#d = d[sample(1:nrow(d), 50000), ]


d$frq = (d$Fk / 5000) * 100

d$est = NA

x = (d$clock == "MutRec")
d$est[x] = d$PostMean[x] * 7300 * (2/100) * 2
x = (d$clock == "Mut")
d$est[x] = d$PostMean[x] * 7300 * (2/100)
x = (d$clock == "Rec")
d$est[x] = d$PostMean[x] * 7300 * 2

x = (d$method == "FGT" | d$method == "FGT phased")
d$est[x] = d$PostMean[x] * 7300 * (1/4)


x = by(d, list(d$clock, d$method), function(x) median(abs(x$time - x$est)) )
array(x, dim(x), dimnames(x))

x = by(d, list(d$clock, d$method), function(x) cor(x$time, x$est, method = "s") )
array(x, dim(x), dimnames(x))



gg = ggplot(d) +
	facet_grid(clock~method) +
	geom_point(aes(x=frq, y=est), colour = "grey50", alpha = 0.05, shape=15, size=1) +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=2.5, colour = "white") +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=1, colour = "grey10") +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=1.5, colour = "blue") +
	coord_cartesian(ylim = c(1, 2000)) +
	scale_y_log10(breaks=c(0, 1, 10, 100, 1000, 10000, 100000)) +
	scale_x_continuous(breaks = (0:5)/10) +
	theme_few() +
	theme(aspect.ratio=2.1/3) +
	xlab("Allele frequency (%)") + ylab("Inferred age (generations)")

ggsave(gg, filename = "_plot.ooa.tru.png", width = 12, height = 10)
ggsave(gg, filename = "_plot.ooa.err.png", width = 12, height = 10)



ggplot(d) +
	facet_grid(clock~method) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, est, colour = Fk), size = 0.5) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.05, 50000), ylim = c(0.05, 50000)) +
	theme(aspect.ratio = 1)  +
	xlab("True age") + ylab("Posterior mode age")





files = dir(pattern = "^t1000\\.truH\\.age\\.pairs\\..+")

p = lapply(files, function(file) {
	x = fread(file, header = T, stringsAsFactors = F)
	x$method = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\1", file)
	x$clock  = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\2", file)
	x
})
p = rbindlist(p)

p = p[which(p$Shared == 1),]
p = p[-(which(p$SegmentLHS == 0)),]
p = p[-(which(p$SegmentRHS == max(p$SegmentRHS))),]

pos = markers$Position
names(pos) = as.character(markers$MarkerID)

p$pos_SegmentLHS = pos[ as.character(p$SegmentLHS) ]
p$pos_SegmentRHS = pos[ as.character(p$SegmentRHS) ]

p$len = (p$pos_SegmentRHS - p$pos_SegmentLHS) + 1

p$str = sprintf("%d %d %d", p$MarkerID + 1, p$SampleID0 + 1, p$SampleID1 + 1)



load("~/Research/DetectIBD/result.truth.local.RData")

truth$str = sprintf("%d %d %d", truth$index, truth$g0, truth$g1)


p = split(p, list(p$method, p$clock))

for (tag in names(p)) {
	sub = p[[tag]]
	
	i = intersect(sub$str, truth$str)
	
	a = match(i, sub$str)
	b = match(i, truth$str)
	
	sub = sub[a, ]
	t = truth[b, ]
	
	sub$len_true = (t$rhs.position - t$lhs.position) + 1
	
	p[[tag]] = sub
}

p = rbindlist(p)



ggplot(p[sample(1:nrow(p), 100000),]) +
	facet_grid(method~clock) +
	geom_point(aes(len_true, len)) +
	scale_x_log10() +
	scale_y_log10()
stat_summary(aes(x=Fk, y=time), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	geom_boxplot(aes(factor(Fk), len))








M = read.table("vanilla.mutations.txt", header=T)
R = read.table("vanilla.records.txt", header=T)

get.age = function(age, M, R) {
	age$node = M$node[ age$MarkerID + 1 ]
	
	x = match(age$node, R$node)
	
	del = which(is.na(x))
	if (length(del) > 0)
		age = age[-del, ]
	
	age$time = R$time[ x ]
	age$length = R$right[ x ] - R$left[ x ]
	
	return(age)
}

read.age = function(file, M, R) {
	age = fread(file, header = T, stringsAsFactors = F)
	get.age(age, M, R)
}






z = fread("x2349.100.age.pairs.FGT.MutRec.HardBreaks.txt", header = T)

z$node = M$node[ z$MarkerID + 1 ]

x = match(z$node[1], R$node)


c0 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\1", R$children))
c1 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\2", R$children))

a = match(60*2+2, c0)
b = match(75*2, c1)
intersect(a, b)


p = data.table(fk  = d$Fk,
							 shr = d$Shared,
							 lhs = M$position[z$SegmentLHS+1],
							 rhs = M$position[z$SegmentRHS+1])
p$len = p$rhs - p$lhs + 1

save(d, file="data.age-ibd.pairs.RData")

ggplot(d) + facet_wrap(~shr) + stat_summary(aes(fk, len), fun.data = "median_hilow", geom = "smooth") + scale_y_log10()


###
stop()
###

library(data.table)

files = dir(pattern="^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$ndiscord = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	tmp$method = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\3", file)
	tmp$clock = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.age.sites.HardBreaks.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$ndiscord = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	tmp$method = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\3", file)
	tmp$clock = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.age.sites.SoftBreaks.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^x([0-9]+)\\.1000\\.age\\.pairs\\.([A-Z]+)\\.MutRec\\.HardBreaks\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$method = sub("^x([0-9]+)\\.([0-9]+)\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\3", file)
	d = rbind(d, tmp)
}

save(d, file="data.age-ibd.pairs.RData")


M = read.table("vanilla.mutations.txt", header=T)

p = data.table(fk  = d$Fk,
							 shr = d$Shared,
							 lhs = M$position[d$SegmentLHS+1],
							 rhs = M$position[d$SegmentRHS+1])
p$len = p$rhs - p$lhs + 1

save(p, file="data.age-ibd.pairs.parsed.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", recursive = T, full.names = T)

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$type   = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
	tmp$method = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\3", file)
	tmp$clock  = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.result.age.sites.RData")




###
ggplot(p) +
	#geom_density_2d(aes(time, PostMode * Ne), bins = 5, size = 1/3, color = "black", alpha = 0.8) +
	#geom_linerange(aes(x = x, ymin = prev, ymax = time)) +
	geom_raster(aes(x, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_raster(aes(mid, parent.time), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_smooth(data = pp, aes(x, Lower * Ne), method = "lm", color = "darkred") +
	#geom_smooth(data = pp, aes(x, Upper * Ne), method = "lm", color = "darkgreen") +
	geom_smooth(data = pp, aes(x, node.time), method = "lm", size = 1/2, color = "grey40") +
	geom_smooth(data = pp, aes(x, mid), method = "lm", size = 1/2, color = "grey40") +
	geom_smooth(data = pp, aes(x, parent.time), method = "lm", size = 1/2, color = "grey40") +
	geom_abline(intercept = c(0,0), slope = 1) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "black", size = 2/3, se = F) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 2/3, fill = "black") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	# geom_point(aes(x, Lower * Ne), alpha=0.25, color = "darkred") +
	# geom_point(aes(x, Upper * Ne), alpha=0.25, color = "darkgreen") +
	# geom_point(aes(x, PostMode * Ne), alpha=0.25) +
	#geom_point(aes(time, prev)) +
	coord_cartesian(xlim = c(0.5, 40000), ylim = c(0.5, 40000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("white", "deepskyblue", "royalblue2", "magenta4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	#scale_fill_gradientn(colours = c("yellow", "orange", "orangered", "darkred", "coral4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("lightcyan2", "deepskyblue", "blue2", "navy"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	#scale_color_gradientn(colours = c("lightblue", "deepskyblue", "royalblue", "darkblue"), na.value = "navy", limits = c(0, 4)) +
	#scale_fill_gradient(low = "green", high = "red") +
	scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				#panel.spacing.y = unit(0.05, "cm"),
				#plot.margin = unit(c(-0.1, 0.05, 0.05, 0.05), "cm"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				#panel.background = element_rect(fill = "grey90"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 8),
				strip.text.y = element_text(face = "bold.italic", angle = 0),
				legend.justification=c(1,0), legend.position=c(0.995,1-0.995), legend.box.background = element_rect(colour = "black"),
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age")
###




###
stop()
###

library(data.table)

files = dir(pattern=".+\\.age\\.sites\\.HMM\\.([a-zA-Z]+)\\.txt$")
files = c(files, dir(pattern=".+\\.age\\.sites\\.HMM\\.([a-zA-Z]+)\\.HardBreaks\\.txt$"))

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	if (grepl("HardBreaks", file)) {
		tmp$type = "Hard breaks"
		tmp$method = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	} else {
		tmp$type = "Soft breaks"
		tmp$method = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	}
	d = rbind(d, tmp)
}

save(d, file="data.result.age.sites.HMM.check.RData")



files = dir(pattern=".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$")
files = c(files, dir(pattern=".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$"))

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	if (grepl("HardBreaks", file)) {
		tmp$type = "Hard breaks"
		tmp$method = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	} else {
		tmp$type = "Soft breaks"
		tmp$method = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	}
	d = rbind(d, tmp)
}

save(d, file="data.result.age.pairs.check.RData")







###############


sim = NULL

d = fread("./truth/simres.age.sites.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Mutation clock"
d$method = "True IBD"
sim = rbind(sim, d)

d = fread("./truth/simres.age.sites.SIM.Rec.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Recombination clock"
d$method = "True IBD"
sim = rbind(sim, d)

d = fread("./truth/simres.age.sites.SIM.MutRec.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Combined clock"
d$method = "True IBD"
sim = rbind(sim, d)


sim = cbind(sim, times[sim$MarkerID + 1, ])


save(sim, file="result.truth.RData")









