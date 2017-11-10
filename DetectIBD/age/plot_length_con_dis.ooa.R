

library(data.table)


marker = fread("./truH.marker.txt", header = T)
times = fread("./OutOfAfricaHapMap20.times.txt", header = T)



TRU = list()
ERR = list()



d = fread("simres.age.pairs.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$tag = "True IBD"
d$method = "True IBD"


d$Fk = d$leaves

d = d[order(d$Fk), ]

d$key = sprintf("%d %d %d %d", d$SampleID0, d$SampleID1, d$SegmentLHS, d$SegmentRHS)

x = lapply(split(d$key, list(d$tag, d$Shared)), duplicated)
sapply(x, table)

d = lapply(split(d, list(d$tag, d$Shared)), function(x) {
	del = which(duplicated(x$key))
	if (length(del) > 0) {
		x = x[-del,]
	}
	x$key = NULL
	x
})
d = rbindlist(d)

table(d$tag, d$Shared)

del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$MarkerID) - 1)
if (length(del) > 0) d = d[-del,]

table(d$tag, d$Shared)


TRU[["TrueIBD"]] = d
ERR[["TrueIBD"]] = d




files = dir(pattern = "^results\\..+_1_0\\.pairs\\.RData$")

for (file in files) {
	print(file)
	
	load(file)
	
	d = cbind(d, times[d$MarkerID + 1, ])
	

	d$Fk = d$leaves
	
	d = d[order(d$Fk), ]
	
	d$key = sprintf("%d %d %d %d", d$SampleID0, d$SampleID1, d$SegmentLHS, d$SegmentRHS)
	
	x = lapply(split(d$key, list(d$tag, d$Shared)), duplicated)
	sapply(x, table)
	
	d = lapply(split(d, list(d$tag, d$Shared)), function(x) {
		del = which(duplicated(x$key))
		if (length(del) > 0) {
			x = x[-del,]
		}
		x$key = NULL
		x
	})
	d = rbindlist(d)
	
	table(d$tag, d$Shared)
	
	del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$MarkerID) - 1)
	if (length(del) > 0) d = d[-del,]
	
	table(d$tag, d$Shared)
	
	
	d$method = ""
	
	if (grepl("fgt_truH", d$tag[1])) d$method = "FGT"
	if (grepl("fgt_truP", d$tag[1])) d$method = "FGT*"
	if (grepl("fgt_errH", d$tag[1])) d$method = "FGT"
	if (grepl("fgt_errP", d$tag[1])) d$method = "FGT*"
	if (grepl("dgt", d$tag[1])) d$method = "DGT"
	if (grepl("hmm", d$tag[1])) d$method = "HMM"
	
	
	if (grepl("tru", file)) TRU[[ file ]] = d
	if (grepl("err", file)) ERR[[ file ]] = d
}



kt = lapply(TRU, function(x) unique(x$MarkerID))
kt = Reduce(intersect, kt)

ke = lapply(ERR, function(x) unique(x$MarkerID))
ke = Reduce(intersect, ke)

key = intersect(kt, ke)

print(length(key))


load("tmp.RData")

TRU = lapply(TRU, function(x) {
	z = which(x$MarkerID %in% key)
	x[z, ]
})

ERR = lapply(ERR, function(x) {
	z = which(x$MarkerID %in% key)
	x[z, ]
})


save(TRU, ERR, file = "_plotdata.lengths.RData")

######


args = commandArgs(T)

load(args[1])


table(d$tag, d$Shared)


d = d[order(d$Fk), ]

d$key = sprintf("%d %d %d %d", d$SampleID0, d$SampleID1, d$SegmentLHS, d$SegmentRHS)

x = lapply(split(d$key, list(d$tag, d$Shared)), duplicated)
sapply(x, table)

d = lapply(split(d, list(d$tag, d$Shared)), function(x) {
	del = which(duplicated(x$key))
	if (length(del) > 0) {
		x = x[-del,]
	}
	x
})
d = rbindlist(d)


table(d$tag, d$Shared)

d$key = NULL


del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$MarkerID) - 1)
if (length(del) > 0) d = d[-del,]

table(d$tag, d$Shared)



d$phy = ((marker$Position[d$SegmentRHS + 1] - marker$Position[d$SegmentLHS + 1]) + 1) / 1e6
d$gen = ((marker$GenDist[d$SegmentRHS + 1] - marker$GenDist[d$SegmentLHS + 1]) + 1e-08) 

x = by(d$phy, list(d$Shared, d$tag), median)
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(d$gen, list(d$Shared, d$tag), median)
array(x, dim = dim(x), dimnames = dimnames(x))




d = split(d, list(d$Fk, d$Shared, d$tag))



tru = lapply(TRU, function(d) {
	d$phy = ((marker$Position[d$SegmentRHS + 1] - marker$Position[d$SegmentLHS + 1]) + 1) / 1e6
	d$gen = ((marker$GenDist[d$SegmentRHS + 1] - marker$GenDist[d$SegmentLHS + 1]) + 1e-08) 
	d = lapply(split(d, list(d$leaves, d$Shared)), function(x) {
		if (nrow(x) > 10000) x[sample(1:nrow(x), 10000), ]
		phy = fivenum(x$phy)
		gen = fivenum(x$gen)
		shr = "Concordant"; if (x$Shared[1] == 0) shr = "Discordant"
		rbind(data.table(fk = x$leaves[1], sh = shr, method = x$method[1], low = phy[2], med = phy[3], upp = phy[4], mode = "A. Physical length"),
					data.table(fk = x$leaves[1], sh = shr, method = x$method[1], low = gen[2], med = gen[3], upp = gen[4], mode = "B. Genetic length"))
	})
	rbindlist(d)
})
tru = rbindlist(tru)


err = lapply(ERR, function(d) {
	d$phy = ((marker$Position[d$SegmentRHS + 1] - marker$Position[d$SegmentLHS + 1]) + 1) / 1e6
	d$gen = ((marker$GenDist[d$SegmentRHS + 1] - marker$GenDist[d$SegmentLHS + 1]) + 1e-08) 
	d = lapply(split(d, list(d$leaves, d$Shared)), function(x) {
		if (nrow(x) > 10000) x[sample(1:nrow(x), 10000), ]
		phy = fivenum(x$phy)
		gen = fivenum(x$gen)
		shr = "Concordant"; if (x$Shared[1] == 0) shr = "Discordant"
		rbind(data.table(fk = x$leaves[1], sh = shr, method = x$method[1], low = phy[2], med = phy[3], upp = phy[4], mode = "A. Physical length"),
					data.table(fk = x$leaves[1], sh = shr, method = x$method[1], low = gen[2], med = gen[3], upp = gen[4], mode = "B. Genetic length"))
	})
	rbindlist(d)
})
err = rbindlist(err)


save(tru, err, file = "_lengths.RData")
save(tru, err, file = "_lengths_reduced.RData")


stop()
#################################



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





d = fread("../test_rvage/simres.age.pairs.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d = cbind(d, times[d$MarkerID + 1, ])
d$tag = "True IBD"


d$Fk = d$leaves

d = d[order(d$Fk), ]

d$key = sprintf("%d %d %d %d", d$SampleID0, d$SampleID1, d$SegmentLHS, d$SegmentRHS)

x = lapply(split(d$key, list(d$tag, d$Shared)), duplicated)
sapply(x, table)

d = lapply(split(d, list(d$tag, d$Shared)), function(x) {
	del = which(duplicated(x$key))
	if (length(del) > 0) {
		x = x[-del,]
	}
	x
})
d = rbindlist(d)


table(d$tag, d$Shared)

d$key = NULL


del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$MarkerID) - 1)
if (length(del) > 0) d = d[-del,]

table(d$tag, d$Shared)




d$phy = ((marker$Position[d$SegmentRHS + 1] - marker$Position[d$SegmentLHS + 1]) + 1) / 1e6
d$gen = ((marker$GenDist[d$SegmentRHS + 1] - marker$GenDist[d$SegmentLHS + 1]) + 1e-08) 

x = by(d$phy, list(d$Shared, d$tag), median)
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(d$gen, list(d$Shared, d$tag), median)
array(x, dim = dim(x), dimnames = dimnames(x))




d = split(d, list(d$Fk, d$Shared, d$tag))

d = lapply(d, function(x) {
	phy = fivenum(x$phy)
	gen = fivenum(x$gen)
	shr = "Concordant"; if (x$Shared[1] == 0) shr = "Discordant"
	rbind(data.table(fk = x$Fk[1], sh = shr, tag = x$tag[1], low = phy[2], med = phy[3], upp = phy[4], mode = "A. Physical length"),
				data.table(fk = x$Fk[1], sh = shr, tag = x$tag[1], low = gen[2], med = gen[3], upp = gen[4], mode = "B. Genetic length"))
})
d = rbindlist(d)


save(d, file = "length.trueIBD_pairs.RData")


stop()
#################################



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


load("_lengths.RData")
load("_lengths_reduced.RData")


p = tru
p = err


p$method = factor(p$method, levels = c("True IBD", "FGT", "FGT*", "DGT", "HMM"), ordered = T)




ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))

len = ggplot(p) + #[sample(1:nrow(p), 10000),]) +
	facet_wrap(~mode, scales = "free", ncol = 1) +
	#geom_violin(aes(fk, x, fill = type), size = 0, alpha = 0.9) +
	#geom_pointrange(data=q, aes(x=fk, y=med, ymin=low, ymax=upp, group=type), color = "black", size = 2/3, shape='|', position=position_dodge(width = 0.9)) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=method, linetype=sh), size = 1.1, position=position_dodge(width = 0.8)) +
	#geom_point(aes(x=fk, y=med, color=type, fill = sh), shape = 18, size = 3, position=position_dodge(width = 0.8)) +
	geom_point(aes(x=fk, y=med, color=method, fill = sh), shape = 23, size = 2, position=position_dodge(width = 0.8)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), minor_breaks = mb, labels = as.character(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10))) +
	scale_x_continuous(breaks = 2:25) +
	scale_color_manual(values = c("grey20", "royalblue1", "purple", "limegreen", "orange3")) +
	#scale_color_manual(values = c("grey20", "royalblue1", "limegreen", "orange3")) +
	scale_linetype_manual(values = rep("solid", 2)) +
	#scale_shape_manual(values = c(18, 20)) +
	scale_fill_manual(values = c("black", "white")) +
	coord_cartesian(ylim = c(0.0005, 29.5), xlim = c(1.5, 25.5), expand = F) +
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

ggsave(len, filename = "___plot.length.con_dis.TRU.pdf", height=6, width = 9)
ggsave(len, filename = "___plot.length.con_dis.ERR.pdf", height=6, width = 9)



