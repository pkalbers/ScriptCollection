
library(data.table)
library(ggplot2)
library(ggthemes)


m = fread("../data/1000G.chr20.marker.txt", header = T, stringsAsFactors = F)


fgt = fread("result.chr20.ibd.FGT.txt", header = T, stringsAsFactors = F)
dgt = fread("result.chr20.ibd.DGT.txt", header = T, stringsAsFactors = F)
hmm = fread("result.chr20.ibd.HMM.txt", header = T, stringsAsFactors = F)


fgt$key = sprintf("%d %d %d %d", fgt$SampleID0, fgt$SampleID1, fgt$LHS, fgt$RHS)
dgt$key = sprintf("%d %d %d %d", dgt$SampleID0, dgt$SampleID1, dgt$LHS, dgt$RHS)
hmm$key = sprintf("%d %d %d %d", hmm$SampleID0, hmm$SampleID1, hmm$LHS, hmm$RHS)


fgt = fgt[order(fgt$Fk, sample(1:nrow(fgt))),]
dgt = dgt[order(dgt$Fk, sample(1:nrow(dgt))),]
hmm = hmm[order(hmm$Fk, sample(1:nrow(hmm))),]


del = which(duplicated(fgt$key));  fgt = fgt[-del,]
del = which(duplicated(dgt$key));  dgt = dgt[-del,]
del = which(duplicated(hmm$key));  hmm = hmm[-del,]



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

p = rbind(fgt, dgt, hmm)

x = sample(1:nrow(fgt), 1e6)
p = rbind(fgt[x,], dgt[x,], hmm[x,])


del = which(p$LHS < 2 | p$RHS > max(p$RHS) - 2)
p = p[-del,]

p$pLHS = m$Position[p$LHS+1]
p$pRHS = m$Position[p$RHS+1]
p$pLen = p$pRHS - p$pLHS + 1

p$gLHS = m$GenDist[p$LHS+1]
p$gRHS = m$GenDist[p$RHS+1]
p$gLen = p$gRHS - p$gLHS + 1e-8


p = split(p, p$tag)
k = Reduce(intersect, sapply(p, function(x) unique(x$key)))
k = sample(k, 1e6)
p = lapply(p, function(x, k) {
	i = which(x$key %in% k)
	x[i,]
}, k)
p = rbindlist(p)



panel = read.table("../data/integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)
panel = as.data.table(panel)


p$s0 = panel$super_pop[p$SampleID0+1]
p$s1 = panel$super_pop[p$SampleID1+1]

p$ss0 = ""
p$ss1 = ""
x = (p$s0 <= p$s1); p$ss0[x] = p$s0[x]; p$ss1[x] = p$s1[x]; 
x = (p$s0 >  p$s1); p$ss0[x] = p$s1[x]; p$ss1[x] = p$s0[x]; 

p$ss = sprintf("%s %s", p$ss0, p$ss1)

nss = by(p, p$ss, nrow)
nss = data.table(ss0 = sub("^([A-Z]+) ([A-Z]+)$", "\\1", names(nss)),
								 ss1 = sub("^([A-Z]+) ([A-Z]+)$", "\\2", names(nss)),
								 num = as.vector(nss))
nss$str = format(nss$num/3, trim = F, big.mark = ",")


p$tag = factor(p$tag, levels = c("fgt", "dgt", "hmm"), labels = c("FGT, phased haplotypes", "DGT, genotypes", "HMM, genotypes"), ordered = T)






gg = ggplot(p) +
	facet_grid(ss1~ss0, switch = "both") +
	stat_summary(aes(x = Fk, y = gLen / 1, color = tag, linetype = tag), geom = "line", size = 2/3, #position = position_dodge(width = 1/2),
							 fun.data = function(x) { data.table(ymin = quantile(x, 1/4),
							 																		y    = quantile(x, 2/4),
							 																		ymax = quantile(x, 3/4)) }) +
	geom_label(data = nss, aes(x = 22, y = 1.65, label = str), label.size = NA, size = 2.5, fill = "white", alpha = 2/3) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10)) +
	scale_colour_manual(values = c("FGT, phased haplotypes" = "purple", "DGT, genotypes" = "limegreen", "HMM, genotypes" = "goldenrod")) +
	scale_linetype_manual(values = c("FGT, phased haplotypes" = "11", "DGT, genotypes" = "22", "HMM, genotypes" = "solid")) +
	coord_cartesian(ylim = c(10^-2.5, 2)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
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
	geom_label(data = nss, aes(x = 22, y = 1.65, label = str), label.size = NA, size = 2.5, fill = "grey90", alpha = 2/3) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10), position = "right") +
	scale_x_continuous(position = "top") +
	scale_colour_manual(values = c("FGT, phased haplotypes" = "purple", "DGT, genotypes" = "limegreen", "HMM, genotypes" = "goldenrod")) +
	scale_linetype_manual(values = c("FGT, phased haplotypes" = "11", "DGT, genotypes" = "22", "HMM, genotypes" = "solid")) +
	coord_cartesian(ylim = c(10^-2.5, 2)) +
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







##################








chr = 20

Ne = 20000


M = fread(sprintf("../../data/1000G.chr%d.marker.txt", chr), header = T, stringsAsFactors = F)


files = dir(pattern = ".*\\.ccf\\..*", full.names = T)

d = NULL

for (file in files) {
	
	if (grepl("NN", file)) next
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	d = rbind(d, tmp)
	
}




loc = unique(d$MarkerID)
length(loc) #  40789
#loc = sample(loc, 5000)


#x = which(d$MarkerID %in% loc)
#p = d[x,]


p = d



x = which(p$Clock == "C")
p = p[x,]
nrow(p)  #  




del = which(p$Shared == 1)
p = p[-del,]
nrow(p)  #  3536825


# key = sprintf("%d %d %d %d %d", p$MarkerID, p$SampleID0, p$Chr0, p$SampleID1, p$Chr1)
# del = which(duplicated(key))
# if (length(del) > 0) {
# 	p = p[-del,]
# }
# nrow(p)  #  1171257
# 
# 
# key = sprintf("%d %d %d %d %d %d", p$SampleID0, p$Chr0, p$SampleID1, p$Chr1, p$SegmentLHS, p$SegmentRHS)
# del = which(duplicated(key))
# if (length(del) > 0) {
# 	p = p[-del,]
# }
# nrow(p)  #  900422


del = which(p$SegmentLHS < 2 | p$SegmentRHS > max(p$SegmentRHS) - 2)
p = p[-del,]
nrow(p)  #  3514601


# x = sample(1:nrow(p), 1e6)
# p = p[x,]
# nrow(p)  #  




p$pLHS = M$Position[p$SegmentLHS+1]
p$pRHS = M$Position[p$SegmentRHS+1]
p$pLen = p$pRHS - p$pLHS + 1

p$gLHS = M$GenDist[p$SegmentLHS+1]
p$gRHS = M$GenDist[p$SegmentRHS+1]
p$gLen = p$gRHS - p$gLHS + 1e-8




panel = read.table("../../data/integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)



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



q = p
key = sprintf("%d %d %d %d %d %d", q$SampleID0, q$Chr0, q$SampleID1, q$Chr1, q$SegmentLHS, q$SegmentRHS)
del = which(duplicated(key))
if (length(del) > 0) {
	q = q[-del,]
}
nrow(q)  #  2338288


x = by(q, list(q$ss), function(x) median(x$gLen))
array(x, dim = dim(x), dimnames = dimnames(x))

y = by(q, list(q$ss), function(x) median(x$pLen/1e6))
array(y, dim = dim(y), dimnames = dimnames(y))

array(x, dim = dim(x), dimnames = dimnames(x)) / array(y, dim = dim(y), dimnames = dimnames(y))




mx = range(round(log10(c(1.925,max(p$Fk)))*(1/3),2))

gg = ggplot(p) +
	facet_grid(ss1~ss0, switch = "both") +
	stat_summary(aes(x = round(log10(Fk)*(1/3),2), y = gLen / 1), geom = "smooth", color = "black", alpha = 1/2, size = 1/2,# fatten = 1/10, 
							 fun.data = function(x) { 
							 	# if (length(x) < 500) {
							 	# 	return(	data.table(ymin = NA,
							 	# 										 y    = NA,
							 	# 										 ymax = NA) )
							 	# }
							 	data.table(ymin = quantile(x, 1/4),
							 						 y    = quantile(x, 2/4),
							 						 ymax = quantile(x, 3/4)) 
							 }) +
	# stat_summary(aes(x = Fk, y = gLen / 1), geom = "pointrange", size = 1/3, fatten = 1/10, color = "indianred",
	# 						 fun.data = function(x) { 
	# 						 	if (length(x) >= 500) {
	# 						 		return(	data.table(ymin = NA,
	# 						 											 y    = NA,
	# 						 											 ymax = NA) )
	# 						 	}
	# 						 	data.table(ymin = quantile(x, 1/4),
	# 						 						 y    = quantile(x, 2/4),
	# 						 						 ymax = quantile(x, 3/4)) 
	# 						 }) +
	geom_label(data = nss, aes(x = 1, y = 2.333, label = str), label.size = NA, size = 2.75, fill = "white", alpha = 2/3) +
	#geom_label(data = avr, aes(x = 41, y = 0.006, label = str), label.size = NA, size = 3, fill = "white", alpha = 2/3, parse = T) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10)) +
	scale_x_continuous(breaks = log10(2^(1:12))*(1/3), labels = c("2","","8","","32","","128","","512","", "2048", ""), expand = c(0,0)) +
	coord_cartesian(ylim = c(0.005, 3), xlim = mx) +
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

ggsave(gg, filename = "___plot.pop_genlen.pdf", width = 12, height = 9)



gg = ggplot(p) +
	facet_grid(ss0~ss1) +
	stat_summary(aes(x = round(log10(Fk)*(1/3),2), y = gLen / 1), geom = "smooth", color = "black", alpha = 1/2, size = 1/2,# fatten = 1/10, 
							 fun.data = function(x) { 
							 	# if (length(x) < 500) {
							 	# 	return(	data.table(ymin = NA,
							 	# 										 y    = NA,
							 	# 										 ymax = NA) )
							 	# }
							 	data.table(ymin = quantile(x, 1/4),
							 						 y    = quantile(x, 2/4),
							 						 ymax = quantile(x, 3/4)) 
							 }) +
	# stat_summary(aes(x = Fk, y = gLen / 1), geom = "pointrange", size = 1/3, fatten = 1/10, color = "indianred",
	# 						 fun.data = function(x) { 
	# 						 	if (length(x) >= 500) {
	# 						 		return(	data.table(ymin = NA,
	# 						 											 y    = NA,
	# 						 											 ymax = NA) )
	# 						 	}
	# 						 	data.table(ymin = quantile(x, 1/4),
	# 						 						 y    = quantile(x, 2/4),
	# 						 						 ymax = quantile(x, 3/4)) 
	# 						 }) +
	#geom_label(data = nss, aes(x = 42, y = 2.33, label = str), label.size = NA, size = 2.75, fill = "grey90", alpha = 2/3) +
	#geom_label(data = mpl, aes(x = 44.5, y = 0.005, label = str), label.size = NA, size = 2.75, fill = "grey90", alpha = 2/3) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10), position = "right") +
	scale_x_continuous(breaks = log10(2^(1:12))*(1/3), labels = c("2","","8","","32","","128","","512","", "2048", ""), expand = c(0,0), position = "top") +
	coord_cartesian(ylim = c(0.005, 3), xlim = mx) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				axis.text = element_text(size = 8),
				strip.placement = "outside",
				strip.text = element_text(face = "bold"),
				panel.background = element_rect(fill = "white"),
				panel.border = element_rect(fill = NA, size = 1/2, colour = "black"),
				panel.grid.major = element_line(colour = "grey90", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey90", size = 1/3)) +
	xlab("Focal allele count") + ylab("Physical length (Mb)")
gg

ggsave(gg, filename = "___plot.pop_phylen.pdf", width = 12, height = 9)




save(p, file = "_data.ibd.hhmm.RData")





##############




















mean(p$gLen / (p$pLen / 1e6))

z = by(p, list(p$ss1, p$ss0), function(x) mean(x$gLen / (x$pLen / 1e6)))
array(z, dim = dim(z), dimnames = dimnames(z))


#############


ggplot(p) +
	stat_summary(aes(x = Fk, y = gLen / 1, color = tag), geom = "pointrange", position = position_dodge(width = 1/2),
							 fun.data = function(x) { data.table(ymin = quantile(x, 1/4),
							 																		y    = quantile(x, 2/4),
							 																		ymax = quantile(x, 3/4)) }) +
	scale_y_log10(minor_breaks = c((0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10))




p$sb = (p$s0 == p$s1)

pp = which(p$sb)
pp = p[pp,]


ggplot(pp) +
	facet_grid(.~s0) +
	stat_summary(aes(x = Fk, y = gLen / 1, color = tag), geom = "pointrange", position = position_dodge(width = 1/2),
							 fun.data = function(x) { data.table(ymin = quantile(x, 1/4),
							 																		y    = quantile(x, 2/4),
							 																		ymax = quantile(x, 3/4)) }) +
	scale_y_log10(minor_breaks = c((0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10))




pp = which(p$s0 != p$s1)
pp = p[pp,]








del = which(duplicated(k))

p = d[-del,]



ggplot(m) +
stat_summary(aes(x = fk, y = len / 1), geom = "pointrange", size = 1/3, fatten = 1/10, #alpha = 1, color = "grey30", #position = position_dodge(width = 1/2),
						 fun.data = function(x) { 
						 	if (length(x) < 100) {
						 		return(	data.table(ymin = NA,
						 											 y    = NA,
						 											 ymax = NA) )
						 	}
						 	data.table(ymin = quantile(x, 1/4),
						 						 y    = quantile(x, 2/4),
						 						 ymax = quantile(x, 3/4)) 
						 }) +
	scale_y_log10(minor_breaks = c((0:9)/1000, (0:9)/100, (0:9)/10, (0:9)/1, (0:9)*10))
