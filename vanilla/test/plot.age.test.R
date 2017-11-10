
library(data.table)
library(ggplot2)
library(ggthemes)


Ne = 1 * 10000

marker = fread("../vanilla.marker.txt", header = T)
times = fread("../vanilla.times.txt", header = T)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}


mx = c(0.5, 20000)
tr = rbind(data.table(method = "   True IBD   ", clock = "(a) Mutation clock",      x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]),
					 data.table(method = "   True IBD   ", clock = "(b) Recombination clock", x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]),
					 data.table(method = "   True IBD   ", clock = "(c) Combined clock",      x0 = mx[1], y0 = mx[1], x1 = mx[2], y1 = mx[2]))

tm = rbind(data.table(x = 10^(0:6), y = 0.001, xend = 10^(0:6), yend = 0.65),
					 data.table(y = 10^(0:6), x = 0.001, yend = 10^(0:6), xend = 0.65))


p = NULL

d = fread("./test.age.sites.SIM.Mut.SoftBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$clock = "(a) Mutation clock"
d$method = "   True IBD   "
p = rbind(p, d)

d = fread("./test.age.sites.SIM.Rec.SoftBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$clock = "(b) Recombination clock"
d$method = "   True IBD   "
p = rbind(p, d)

d = fread("./test.age.sites.SIM.MutRec.SoftBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$clock = "(c) Combined clock"
d$method = "   True IBD   "
p = rbind(p, d)


p$mid = exp(log(p$node.time) + ((log(p$parent.time) - log(p$node.time)) / 2))
p$x = p$mid # p$node.time

pp = p[which(p$node.time >= 1 & p$node.time <= 1000), ]


ggplot(p) +
	facet_grid(method~clock) +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	geom_raster(aes(x, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(10/100, 10/100)) +
	geom_smooth(data = pp, aes(x, node.time), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	#geom_smooth(data = pp, aes(x, mid), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	geom_smooth(data = pp, aes(x, parent.time), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	geom_abline(intercept = c(0,0), slope = 1) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = T) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", fullrange = T) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_fill_gradientn(colours = c("lightcyan2", "lightblue2", "deepskyblue1", "royalblue1", "purple1", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 9),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				#legend.justification=c(1,0), legend.position=c(0.995,1-0.985), legend.box.background = element_rect(colour = "black"),
				legend.position = "bottom",
				#legend.key.height = unit(0.4, "cm"),
				#legend.key.width = unit(0.3, "cm"),
				legend.key.width = unit(1, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age")



z = sapply(split(p, list(p$clock, p$method)), function(x) {
	data.table("Child node"  = cor(x$PostMode * Ne, x$node.time, method = "s"),
						 "Midpoint"    = cor(x$PostMode * Ne, x$mid, method = "s"),
						 "Parent node" = cor(x$PostMode * Ne, x$parent.time, method = "s"))
	
})
t(array(z, dim = dim(z), dimnames = dimnames(z)))


z = sapply(split(p, list(p$clock, p$method)), function(x) {
	data.table("Child node"  = rmsle(x$PostMode * Ne, x$node.time),
						 "Midpoint"    = rmsle(x$PostMode * Ne, x$mid),
						 "Parent node" = rmsle(x$PostMode * Ne, x$parent.time))
	
})
t(array(z, dim = dim(z), dimnames = dimnames(z)))





a = fread("./test.age.sites.SIM.Rec.SoftBreaks.txt", header = T, stringsAsFactors = F)
b = fread("../simres2.age.sites.SIM.Rec.HardBreaks.txt", header = T, stringsAsFactors = F)
i = intersect(a$MarkerID, b$MarkerID)
a = a[which(a$MarkerID %in% i),]
b = b[which(b$MarkerID %in% i),]
a = a[order(a$MarkerID),]
b = b[order(b$MarkerID),]

plot(density(b$PostMode))
lines(density(a$PostMode), col = "red")
