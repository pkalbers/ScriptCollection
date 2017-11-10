

library(data.table)
library(ggplot2)
library(ggthemes)


d = list()

for (i in 1:22) {
	print(i)
	
	M = fread(sprintf("../data/1000G.chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	pl = (max(M$Position) - min(M$Position)) + 1
	gl = (max(M$GenDist) - min(M$GenDist))
	
	fgt = read.table(sprintf("./length.chr%d.ibd.FGT.txt", i), header = T, stringsAsFactors = F)
	dhg = read.table(sprintf("./length.chr%d.ibd.DHG.txt", i), header = T, stringsAsFactors = F)
	
	fgt$method = " Four-gamete test "
	dhg$method = "Discordant homozygote genotypes"
	
	key = sprintf("%d %d", fgt$SampleID0, fgt$SampleID1); del = which(duplicated(key)); if (length(del) > 0) fgt = fgt[-del,]
	key = sprintf("%d %d", dhg$SampleID0, dhg$SampleID1); del = which(duplicated(key)); if (length(del) > 0) dhg = dhg[-del,]
	
	tmp = rbind(fgt, dhg)
	
	tmp$pos.lhs = M$Position[ tmp$LHS + 1 ]
	tmp$pos.rhs = M$Position[ tmp$RHS + 1 ]
	tmp$pos.len = (tmp$pos.rhs - tmp$pos.lhs) + 1
	tmp$pos.frac = tmp$pos.len / pl
	
	tmp$gen.lhs = M$GenDist[ tmp$LHS + 1 ]
	tmp$gen.rhs = M$GenDist[ tmp$RHS + 1 ]
	tmp$gen.len = (tmp$gen.rhs - tmp$gen.lhs)
	tmp$gen.frac = tmp$gen.len / gl
	
	tmp$chr = i
	
	tmp$wall = F
	tmp$wall[which(tmp$LHS < 100 | tmp$RHS > max(M$MarkerID) - 100)] = T
	
	d[[as.character(i)]] = tmp
}


sapply(d, nrow)
sapply(d, function(x) table(x$Fk))
sapply(d, function(x) table(x$method))


d = lapply(d, function(d) {
	del = which(d$Fk == 26)
	if (length(del) > 0)
		d = d[-del,]
	
	del = which(is.na(d$pos.len))
	if (length(del) > 0)
		d = d[-del,]
	
	del = which(is.na(d$gen.len))
	if (length(del) > 0)
		d = d[-del,]
	
	del = which(d$wall)
	if (length(del) > 0)
		d = d[-del,]
	
	d
})


sapply(d, nrow)
sapply(d, function(x) table(x$Fk))
sapply(d, function(x) table(x$method))



d = rbindlist(d)

table(d$method)


a = sort(unique(d$Fk))
b = (a / 5000) * 100
frq = b
names(frq) = a

x = ((a %% 5 == 0))
frq[ !x ] = ""


gg = ggplot(d) + 
	facet_grid(.~method) + 
	geom_boxplot(aes(x = factor(Fk), y = pos.len), fill = "grey70", outlier.shape = NA) + #outlier.colour = "white", outlier.size = 0) +
	scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 3e7, by = 0.02e7), labels = sprintf("%.1f", seq(0, 3e7, by = 0.02e7) / 1e6)) +
	scale_x_discrete(labels = frq) +
	coord_cartesian(ylim=c(0, 0.2e7)) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold"),
				panel.grid = element_line(size = 0.5, colour = "grey"),
				panel.grid.major = element_line(size = 0.5, colour = "grey"),
				panel.grid.major.y = element_line(size = 0.5, colour = "grey"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank(),
				panel.grid.minor.x = element_blank()) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (Mb)")
gg
ggsave(gg, filename = "_plot.phylength.pdf", height = 5, width = 10)
ggsave(gg, filename = "_plot.phylength.png", height = 5, width = 10)


gg = ggplot(d) + 
	facet_grid(.~method) + 
	geom_boxplot(aes(x = factor(Fk), y = gen.len), fill = "grey70", outlier.shape = NA) + #outlier.colour = "white", outlier.size = 0) +
	scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 3, by = 0.2), labels = sprintf("%.1f", seq(0, 3, by = 0.2))) +
	scale_x_discrete(labels = frq) +
	coord_cartesian(ylim=c(0, 2)) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (cM)")
gg
ggsave(gg, filename = "_plot.genlength.pdf", height = 5, width = 10)
ggsave(gg, filename = "_plot.genlength.png", height = 5, width = 10)



by(d, d$method, function(x) median(x$pos.len, na.rm = T))
by(d, d$method, function(x) median(x$gen.len, na.rm = T))




col = c("dodgerblue2","#E31A1C", # red
	"green4",
	"#6A3D9A", # purple
	"#FF7F00", # orange
	"black", #"gold1",
	"skyblue2","#FB9A99", # lt pink
	"palegreen2",
	"#CAB2D6", # lt purple
	"#FDBF6F", # lt orange
	"gray70", #"khaki2",
	"maroon","orchid1","deeppink1","blue1","steelblue4",
	"darkturquoise","green1","yellow4","yellow3",
	"darkorange4","brown")

ggplot(d) + 
	facet_grid(.~method) + 
	stat_summary(aes(x = (Fk), y = gen.len, colour = factor(chr)), fun.y = "median", geom="line") + 
	coord_cartesian(ylim = c(0, 0.55)) + 
	scale_colour_manual(values = col) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (cM)")


x = by(d, list(d$method, d$chr), function(x) round(median(x$pos.len / 1e6, na.rm = T), 5))
x = t(array(x, dim(x), dimnames = dimnames(x)))

y = by(d, list(d$method, d$chr), function(x) round(median(x$gen.len, na.rm = T), 5))
y = t(array(y, dim(y), dimnames = dimnames(y)))

r = by(d, list(d$method, d$chr), function(x) round(mean(x$gen.len / (x$pos.len/1e6), na.rm = T), 5))
r = t(array(r, dim(r), dimnames = dimnames(r)))

se <- function(x) sqrt(var(x)/length(x))

s = by(d, list(d$method, d$chr), function(x) round(se(x$gen.len / (x$pos.len/1e6)), 5))
s = t(array(s, dim(s), dimnames = dimnames(s)))


a = sprintf("%.2f (%.3f)", r[,1], s[,1])
b = sprintf("%.2f (%.3f)", r[,2], s[,2])


x = apply(x, 2, function(x) sprintf("%.2f", x))
y = apply(y, 2, function(x) sprintf("%.2f", x))

z = cbind(y, x, a, b)
z = cbind(1:22, z)
colnames(z) = c("Chr.", "FGT", "DHG", "FGT", "DHG", "FGT", "DHG")

write.csv(z, file = "_table.length_chr.csv", quote = F, row.names = T)







