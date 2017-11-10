
library(data.table)
library(ggplot2)
library(ggthemes)


m = fread("../data/1000G.chr20.marker.txt", header = T, stringsAsFactors = F)


fgt = fread("result.chr20.ibd.FGT.txt", header = T, stringsAsFactors = F)
dgt = fread("result.chr20.ibd.DGT.txt", header = T, stringsAsFactors = F)
hmm = fread("result.chr20.ibd.HMM.txt", header = T, stringsAsFactors = F)


# fgt$key = sprintf("%d %d %d %d", fgt$SampleID0, fgt$SampleID1, fgt$LHS, fgt$RHS)
# dgt$key = sprintf("%d %d %d %d", dgt$SampleID0, dgt$SampleID1, dgt$LHS, dgt$RHS)
# hmm$key = sprintf("%d %d %d %d", hmm$SampleID0, hmm$SampleID1, hmm$LHS, hmm$RHS)
# 
# 
# fgt = fgt[order(fgt$Fk, sample(1:nrow(fgt))),]
# dgt = dgt[order(dgt$Fk, sample(1:nrow(dgt))),]
# hmm = hmm[order(hmm$Fk, sample(1:nrow(hmm))),]
# 
# 
# del = which(duplicated(fgt$key));  fgt = fgt[-del,]
# del = which(duplicated(dgt$key));  dgt = dgt[-del,]
# del = which(duplicated(hmm$key));  hmm = hmm[-del,]
# 


rhs = max(c(fgt$RHS, dgt$RHS, hmm$RHS))

del = which(fgt$LHS < 2 | fgt$RHS > rhs - 2); if (length(del) > 0) fgt = fgt[-del,]
del = which(dgt$LHS < 2 | dgt$RHS > rhs - 2); if (length(del) > 0) dgt = dgt[-del,]
del = which(hmm$LHS < 2 | hmm$RHS > rhs - 2); if (length(del) > 0) hmm = hmm[-del,]



fgt$key = sprintf("%d %d %d", fgt$MarkerID, fgt$SampleID0, fgt$SampleID1)
dgt$key = sprintf("%d %d %d", dgt$MarkerID, dgt$SampleID0, dgt$SampleID1)
hmm$key = sprintf("%d %d %d", hmm$MarkerID, hmm$SampleID0, hmm$SampleID1)

key = intersect(hmm$key, intersect(fgt$key, dgt$key))

fgt = as.data.table(as.data.frame(fgt)[match(key, fgt$key),])
dgt = as.data.table(as.data.frame(dgt)[match(key, dgt$key),])
hmm = as.data.table(as.data.frame(hmm)[match(key, hmm$key),])


fgt$tag = "fgt"
dgt$tag = "dgt"
hmm$tag = "hmm"

#p = rbind(fgt, dgt, hmm)


fgt$key = sprintf("%d %d %d %d", fgt$SampleID0, fgt$SampleID1, fgt$LHS, fgt$RHS)
dgt$key = sprintf("%d %d %d %d", dgt$SampleID0, dgt$SampleID1, dgt$LHS, dgt$RHS)
hmm$key = sprintf("%d %d %d %d", hmm$SampleID0, hmm$SampleID1, hmm$LHS, hmm$RHS)

fgt = fgt[order(fgt$MarkerID, fgt$SampleID0, fgt$SampleID1),]
dgt = dgt[order(dgt$MarkerID, dgt$SampleID0, dgt$SampleID1),]
hmm = hmm[order(hmm$MarkerID, hmm$SampleID0, hmm$SampleID1),]

identical(fgt$MarkerID, dgt$MarkerID)
identical(fgt$MarkerID, hmm$MarkerID)
identical(fgt$SampleID0, dgt$SampleID0)

x = sample(1:nrow(fgt), 2e6)
p = rbind(fgt[x,], dgt[x,], hmm[x,])




p$pLHS = m$Position[p$LHS+1]
p$pRHS = m$Position[p$RHS+1]
p$pLen = p$pRHS - p$pLHS + 1

p$gLHS = m$GenDist[p$LHS+1]
p$gRHS = m$GenDist[p$RHS+1]
p$gLen = p$gRHS - p$gLHS + 1e-8


# p = split(p, p$tag)
# k = Reduce(intersect, sapply(p, function(x) unique(x$key)))
# k = sample(k, 1e6)
# p = lapply(p, function(x, k) {
# 	i = which(x$key %in% k)
# 	x[i,]
# }, k)
# p = rbindlist(p)



panel = read.table("../data/integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)
panel = as.data.table(panel)


p$s0 = panel$super_pop[p$SampleID0+1]
p$s1 = panel$super_pop[p$SampleID1+1]

p$ss0 = ""
p$ss1 = ""
x = (p$s0 <= p$s1); p$ss0[x] = p$s0[x]; p$ss1[x] = p$s1[x]; 
x = (p$s0 >  p$s1); p$ss0[x] = p$s1[x]; p$ss1[x] = p$s0[x]; 

p$ss = sprintf("%s %s", p$ss0, p$ss1)



nss = table(p$ss)
nss = data.table(ss0 = sub("^([A-Z]+) ([A-Z]+)$", "\\1", names(nss)),
								 ss1 = sub("^([A-Z]+) ([A-Z]+)$", "\\2", names(nss)),
								 num = as.vector(nss))
nss$str = sprintf("%s %%", format(round(nss$num/sum(nss$num) * 100, 3), digits = 3, trim = F, big.mark = ","))



q = split(p, p$ss)
q = lapply(q, function(x) {
	x = x[order(x$Fk, sample(1:nrow(x))),]
	
	x = split(x, x$tag)

	x = lapply(x, function(y) {
		del = which(duplicated(y$key))
		if (length(del) > 0) {
			y = y[-del,]
		}
		y
	})
	
	k = Reduce(intersect, lapply(x, function(y) { unique(y$MarkerID) }))
	
	x = lapply(x, function(y, k) {
		i = which(y$MarkerID %in% k)
		y[i,]
	}, k)
	
	rbindlist(x)
})
p = rbindlist(q)



p$tag = factor(p$tag, levels = c("fgt", "dgt", "hmm"), labels = c("FGT, phased haplotypes", "DGT, genotypes", "HMM, genotypes"), ordered = T)






gg = ggplot(p) +
	facet_grid(ss1~ss0, switch = "both") +
	stat_summary(aes(x = Fk, y = gLen / 1, color = tag, linetype = tag), geom = "line", size = 2/3, #position = position_dodge(width = 1/2),
							 fun.data = function(x) { data.table(ymin = quantile(x, 1/4),
							 																		y    = quantile(x, 2/4),
							 																		ymax = quantile(x, 3/4)) }) +
	geom_label(data = nss, aes(x = 21.5, y = 2.333, label = str), label.size = NA, size = 2.5, fill = "white", alpha = 2/3) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10)) +
	scale_colour_manual(values = c("FGT, phased haplotypes" = "purple", "DGT, genotypes" = "limegreen", "HMM, genotypes" = "goldenrod")) +
	scale_linetype_manual(values = c("FGT, phased haplotypes" = "11", "DGT, genotypes" = "22", "HMM, genotypes" = "solid")) +
	coord_cartesian(ylim = c(0.005, 3)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				axis.text = element_text(size = 8),
				strip.placement = "outside",
				strip.text = element_text(face = "bold"),
				panel.background = element_rect(fill = "grey90"),
				panel.border = element_rect(fill = NA, size = 1/2, colour = "black"),
				panel.grid.major = element_line(colour = "white", size = 1/2),
				panel.grid.minor.y = element_line(colour = "white", size = 1/3)) +
	xlab("Focal allele count") + ylab("Genetic length (cM)")
gg

ggsave(gg, filename = "__plot.pop_genlen.pdf", width = 12, height = 9)



gg = ggplot(p) +
	facet_grid(ss0~ss1) +
	stat_summary(aes(x = Fk, y = pLen / 1e6, color = tag, linetype = tag), geom = "line", size = 2/3, #position = position_dodge(width = 1/2),
							 fun.data = function(x) { data.table(ymin = quantile(x, 1/4),
							 																		y    = quantile(x, 2/4),
							 																		ymax = quantile(x, 3/4)) }) +
	geom_label(data = nss, aes(x = 21.5, y = 2.333, label = str), label.size = NA, size = 2.5, fill = "grey90", alpha = 2/3) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10), position = "right") +
	scale_x_continuous(position = "top") +
	scale_colour_manual(values = c("FGT, phased haplotypes" = "purple", "DGT, genotypes" = "limegreen", "HMM, genotypes" = "goldenrod")) +
	scale_linetype_manual(values = c("FGT, phased haplotypes" = "11", "DGT, genotypes" = "22", "HMM, genotypes" = "solid")) +
	coord_cartesian(ylim = c(0.005, 3)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				strip.placement = "outside",
				strip.text = element_text(face = "bold"),
				panel.background = element_rect(fill = "white"),
				panel.border = element_rect(fill = NA, size = 1/2, colour = "black"),
				panel.grid.major = element_line(colour = "grey90", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey90", size = 1/3)) +
	xlab("Focal allele count") + ylab("Physical length (Mb)")
gg

ggsave(gg, filename = "__plot.pop_phylen.pdf", width = 12, height = 9)




save(p, file = "_data.ibd.fgt-dgt-hmm.RData")
load("_data.ibd.fgt-dgt-hmm.RData")







x = by(p, list(p$ss, p$tag), function(x) median(x$gLen))
array(x, dim = dim(x), dimnames = dimnames(x))

y = by(p, list(p$ss, p$tag), function(x) median(x$pLen/1e6))
array(y, dim = dim(y), dimnames = dimnames(y))



#############



load("_data.ibd.fgt-dgt-hmm.RData")


q = split(p, p$tag)

q = rbind(data.table(idx = q$`FGT, phased haplotypes`$MarkerID, fk = q$`FGT, phased haplotypes`$Fk, type = "FGT, phased haplotypes", x = q$`FGT, phased haplotypes`$pLen/1e6, mode = "A. Physical length"),
					data.table(idx = q$`FGT, phased haplotypes`$MarkerID, fk = q$`FGT, phased haplotypes`$Fk, type = "FGT, phased haplotypes", x = q$`FGT, phased haplotypes`$gLen, mode = "B. Genetic length"),
					data.table(idx = q$`DGT, genotypes`$MarkerID, fk = q$`DGT, genotypes`$Fk, type = "DGT, genotypes", x = q$`DGT, genotypes`$pLen/1e6, mode = "A. Physical length"),
					data.table(idx = q$`DGT, genotypes`$MarkerID, fk = q$`DGT, genotypes`$Fk, type = "DGT, genotypes", x = q$`DGT, genotypes`$gLen, mode = "B. Genetic length"),
					data.table(idx = q$`HMM, genotypes`$MarkerID, fk = p$Fk, type = "HMM, genotypes", x = q$`HMM, genotypes`$pLen/1e6, mode = "A. Physical length"),
					data.table(idx = q$`HMM, genotypes`$MarkerID, fk = p$Fk, type = "HMM, genotypes", x = q$`HMM, genotypes`$gLen, mode = "B. Genetic length"))

q$type = factor(q$type, levels = c("FGT, phased haplotypes", "DGT, genotypes", "HMM, genotypes"), ordered = T)




x = by(q$x, list(q$type, q$mode), median)
array(x, dim(x), dimnames(x))

x = by(q$x, list(q$fk, q$type, q$mode), median)
array(x, dim(x), dimnames(x))


q = split(q, list(q$fk, q$type, q$mode))
q = lapply(q, function(x) {
	z = x
	w = quantile(x$x, c(1/4,2/4,3/4))
	z$low = w[1]
	z$med = w[2]
	z$upp = w[3]
	z[1,]
})
q = rbindlist(q)


tl = c("FGT, phased haplotypes", "DGT, genotypes", "HMM, genotypes")

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
	scale_color_manual(values = c("purple", "limegreen", "goldenrod")) +
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

ggsave(len, filename = "___plot.chr20.length.fgt-dgt-hmm.pdf", height=7, width = 9)







