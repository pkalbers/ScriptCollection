

library(data.table)
library(ggplot2)
library(ggthemes)


Ne = 2 * 10000

marker = fread("./vanilla.marker.txt", header = T)
times = fread("./vanilla.times.txt", header = T)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}






p = NULL

d = fread("./simres.age.pairs.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$method = "True IBD"
p = rbind(p, d)



d = fread("./detect.age.pairs.FGT.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$method = "FGT"
p = rbind(p, d)

d = fread("./detect.age.pairs.DGT.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$method = "DGT"
p = rbind(p, d)

d = fread("./detect.age.pairs.HMM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$method = "HMM"
p = rbind(p, d)



p$method = factor(p$method, levels = c("True IBD", "FGT", "DGT", "HMM"), ordered = T)



save(p, file = "results.discord_length.RData")
load("results.discord_length.RData")





table(p$method, p$Shared)


p$Fk = p$leaves

p = p[order(p$Fk), ]

p$key = sprintf("%d %d %d %d", p$SampleID0, p$SampleID1, p$SegmentLHS, p$SegmentRHS)

x = lapply(split(p$key, list(p$method, p$Shared)), duplicated)
sapply(x, table)

p = lapply(split(p, list(p$method, p$Shared)), function(x) {
	del = which(duplicated(x$key))
	if (length(del) > 0) {
		x = x[-del,]
	}
	x
})
p = rbindlist(p)


table(p$Shared, p$method)



del = which(p$SegmentLHS <= 1 | p$SegmentRHS >= max(marker$MarkerID) - 1)
if (length(del) > 0) p = p[-del,]

table(p$Shared, p$method)



key = split(p$MarkerID, list(p$method, p$Shared))
key = lapply(key, unique)
key = Reduce(intersect, key)

p = split(p, list(p$method, p$Shared))
p = lapply(p, function(x) {
	z = which(x$MarkerID %in% key)
	x[z, ]
})
p = rbindlist(p)

table(p$Shared, p$method)



p$phy = ((marker$Position[p$SegmentRHS + 1] - marker$Position[p$SegmentLHS + 1]) + 1) / 1e6
p$gen = ((marker$GenDist[p$SegmentRHS + 1] - marker$GenDist[p$SegmentLHS + 1]) + 1e-08) 

x = by(p$phy, list(p$Shared, p$method), median)
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(p$gen, list(p$Shared, p$method), median)
array(x, dim = dim(x), dimnames = dimnames(x))




q = split(p, list(p$Fk, p$Shared, p$method))

q = lapply(q, function(x) {
	phy = fivenum(x$phy)
	gen = fivenum(x$gen)
	shr = "Concordant"; if (x$Shared[1] == 0) shr = "Discordant"
	rbind(data.table(fk = x$Fk[1], sh = shr, type = x$method[1], low = phy[2], med = phy[3], upp = phy[4], mode = "A. Physical length"),
				data.table(fk = x$Fk[1], sh = shr, type = x$method[1], low = gen[2], med = gen[3], upp = gen[4], mode = "B. Genetic length"))
})
q = rbindlist(q)




tl = c("True IBD", "FGT", "DGT", "HMM")

q$type = factor(q$type, levels = tl, ordered = T)

ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))

len = ggplot(q) + #[sample(1:nrow(p), 10000),]) +
	facet_wrap(~mode, scales = "free", ncol = 1) +
	#geom_violin(aes(fk, x, fill = type), size = 0, alpha = 0.9) +
	#geom_pointrange(data=q, aes(x=fk, y=med, ymin=low, ymax=upp, group=type), color = "black", size = 2/3, shape='|', position=position_dodge(width = 0.9)) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=type, linetype=sh), size = 1.1, position=position_dodge(width = 0.8)) +
	#geom_point(aes(x=fk, y=med, color=type, fill = sh), shape = 18, size = 3, position=position_dodge(width = 0.8)) +
	geom_point(aes(x=fk, y=med, color=type, fill = sh), shape = 23, size = 2, position=position_dodge(width = 0.8)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), minor_breaks = mb, labels = as.character(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10))) +
	scale_x_continuous(breaks = 2:25) +
	#scale_color_manual(values = c("grey20", "royalblue1", "purple", "limegreen", "orange3")) +
	scale_color_manual(values = c("grey20", "royalblue1", "limegreen", "orange3")) +
	scale_linetype_manual(values = rep("solid", 2)) +
	#scale_shape_manual(values = c(18, 20)) +
	scale_fill_manual(values = c("black", "white")) +
	coord_cartesian(ylim = c(0.0035, 11.5), xlim = c(1.5, 20.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "horizontal",
				legend.margin = margin(1,2,1,0, "mm"), 
				legend.box.margin = margin(0,0,0,0, "mm"), legend.box = "horizontal",
				legend.text = element_text(size = 9),
				legend.key.height = unit(0.5, "cm"),
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
	ylab("Genetic length (cM)                                Physical length (Mb)") +
	guides(color = guide_legend(override.aes = list(size = 3, shape = 32)))
len

ggsave(len, filename = "___plot.length.con_dis.pdf", height=6, width = 9)



