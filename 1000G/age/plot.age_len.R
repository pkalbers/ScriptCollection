

library(ggplot2)
library(ggthemes)
library(data.table)


marker = fread("../data/1000G.chr20.marker.txt", header = T)


file = "result.age.pairs.HMM.Rec.HardBreaks.txt"


tmp = fread(file, header = T, stringsAsFactors = F)

tmp$tmrca = sub(".+_t([0-9]+).+", "\\1", file)

tmp$i.len = (tmp$SegmentRHS - tmp$SegmentLHS) + 1
tmp$p.len = (marker$Position[tmp$SegmentRHS+1] - marker$Position[tmp$SegmentLHS+1]) + 1
tmp$g.len = (marker$GenDist[tmp$SegmentRHS+1] - marker$GenDist[tmp$SegmentLHS+1]) + 1e-8

d = tmp



# for (file in dir(pattern = "^res_.+\\.age\\.pairs\\..+txt$", path = "./res", full.names = T)) {
# 	
# 	tmp = fread(file, header = T, stringsAsFactors = F)
# 	
# 	tmp = cbind(tmp, times[tmp$MarkerID + 1, ])
# 	
# 	tmp$tmrca = sub(".+_t([0-9]+).+", "\\1", file)
# 	tmp$error = F
# 	tmp$phased = "Correct phase"
# 	tmp$clock = ""
# 	tmp$method = ""
# 	
# 	if (grepl("_err", file)) tmp$error = T
# 	
# 	if (grepl("_errP", file) || grepl("_truP", file)) tmp$phased = "Phased"
# 	
# 	if (grepl(".Mut.", file)) tmp$clock = "(a) Mutation clock"
# 	if (grepl(".Rec.", file)) tmp$clock = "(b) Recombination clock"
# 	if (grepl(".MutRec.", file)) tmp$clock = "(c) Combined clock"
# 	
# 	if (grepl(".HMM.", file)) tmp$method = "H-HMM"
# 	if (grepl(".HMMi.", file)) tmp$method = "H-HMM*"
# 	
# 	tmp$i.len = (tmp$SegmentRHS - tmp$SegmentLHS) + 1
# 	tmp$p.len = (marker$Position[tmp$SegmentRHS+1] - marker$Position[tmp$SegmentLHS+1]) + 1
# 	tmp$g.len = (marker$GenDist[tmp$SegmentRHS+1] - marker$GenDist[tmp$SegmentLHS+1]) + 1e-8
# 	
# 	d = rbind(d, tmp)
# }


d$p.len = d$p.len / 1e6

d$cordant = "Discordant"
d$cordant[which(d$Shared == 1)] = "Concordant"


del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$SegmentRHS)-1)
if (length(del) > 0) d = d[-del,]




p = d
p = p[which(!p$error),]
p = p[which(p$method == "H-HMM"),]
p = p[which(p$tmrca == "10000"),]


gr = NULL
for (k in 2:51) gr = rbind(gr, data.table(x=k-0.5, y = 10^(-8:5)))


len = ggplot(p) +
	facet_grid(.~cordant) +
	geom_raster(aes(Fk,    p.len), position = "identity", stat = "bin2d", binwidth = c(1/50,1/10)) +
	#geom_tile(aes(Fk,      g.len), position = "identity", stat = "bin2d", binwidth = c(1/25,1/10)) +
	geom_vline(xintercept = (3:50)-0.5, alpha = 1/10, color = "black") +
	geom_point(data = gr, aes(x, y), alpha = 1/2, color = "black", shape = 95, size = 2) +
	geom_linerange(aes(Fk, p.len), stat = "summary", size = 2, alpha = 1/2, color = "black", fun.ymin = function(z) {quantile(z,0.25)}, fun.ymax = function(z) {quantile(z,0.75)}) +
	geom_point(aes(Fk,     p.len), stat = "summary", alpha = 2/3, color = "black", fun.y = median, shape=19) +
	geom_point(aes(Fk,     p.len), stat = "summary", alpha = 2/3, color = "white", fun.y = median, shape=20) +
	coord_cartesian(xlim = c(1.5, 50.5), ylim = c(7e-09, 150), expand = F) +
	scale_x_continuous(breaks = c(2, 5, 10, 20, 30, 40, 50)) +
	scale_y_log10(breaks = 10^(-8:5)) + # , labels = trimws(format(10^(-8:5), scientific = T))) +
	scale_fill_gradientn(colours = c("grey90", "aquamarine", "aquamarine2", "goldenrod1", "orange", "red", "darkred"), na.value = "darkred", trans = "log10", limits=c(1, 1e4)) +
	theme_few() +
	theme(aspect.ratio = 0.5,
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "grey50", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,0), legend.position=c(1-0.985,1-0.985),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2),
				#legend.position = "bottom",
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	#ylab("Genetic length (cM)") + 
	ylab("Physical length (Mb)") + 
	xlab("Allele count, k")
len

ggsave(len, filename = sprintf("_plot.len-g.hhmm_.t%s.pdf", "100"), height = 3.5, width = 11)
ggsave(len, filename = sprintf("_plot.len-p.hhmm_.t%s.pdf", "100"), height = 3.5, width = 11)
#ggsave(len, filename = "_plot.len-g.hhmmi.pdf", height = 7, width = 5.5)
#ggsave(len, filename = "_plot.len-p.hhmmi.pdf", height = 7, width = 5.5)






### age by frq


file = "result.age.sites.HMM.Rec.HardBreaks.txt"

tmp = fread(file, header = T, stringsAsFactors = F)

tmp$clock = ""
if (grepl(".Mut.", file)) tmp$clock = "(a) Mutation clock"
if (grepl(".Rec.", file)) tmp$clock = "(b) Recombination clock"
if (grepl(".MutRec.", file)) tmp$clock = "(c) Combined clock"

d = tmp



p = d


p$x = p$Fk 



my = c(0.9, 1.2e5)
mx = c(2, 51)

tm = rbind(data.table(y = 10^(0:6), x = 0, yend = 10^(0:6), xend = 1.5))


scat = ggplot(p) +
	#facet_grid(.~clock) +
	#geom_raster(aes(x, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(1/50, 1/5)) +
	geom_tile(aes(x+1, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(1, 1/5)) +
	geom_vline(xintercept = (2:50), color = "black", alpha = 1/10, size = 1/2) +
	#geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = T) +
	#geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", fullrange = T) +
	#geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	coord_cartesian(xlim = mx, ylim = my, expand = F) +
	scale_x_continuous(breaks = c(2, 5, 10, 20, 30, 40, 50)+0.5, labels=c(2, 5, 10, 20, 30, 40, 50)) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy", "black"), na.value = "black", trans = "log10", limits = c(1, 50)) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1/2,
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "grey50", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,1), legend.position=c(1-0.99,0.98),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
				#legend.position = "bottom",
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("Allele count, k")
scat


ggsave(scat, filename = sprintf("_plot.scat.hhmm_.t%s.pdf", "100"), height = 5, width = 10)





