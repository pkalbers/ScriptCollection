

library(data.table)
library(ggplot2)
library(ggthemes)



unique.segments = T
equal.nsegments = F
adjust.breakpts = T
subsample.lengths = T



data = NULL

files = dir(pattern = "length\\.chr[0-9]+\\.ibd\\..+\\.txt$")

for (file in files) {
	print(file)
	tmp = read.table(file, header = T, stringsAsFactors = F)
	
	chr = sub("^length\\.chr([0-9]+)\\.ibd\\..+\\.txt$", "\\1", file)
	
	tmp$chr = chr
	
	
	if (sub("^length\\.chr[0-9]+\\.ibd\\.(.+)\\.txt$", "\\1", file) == "FGT") {
		tmp$type = "(a) FGT, phased haplotypes"
	}
	
	if (sub("^length\\.chr[0-9]+\\.ibd\\.(.+)\\.txt$", "\\1", file) == "DHG") {
		tmp$type = "(b) DGT, genotypes"
	}
	
	
	mrk = sprintf("../data/1000G.chr%s_c.marker.txt", chr)
	mrk = fread(mrk, header = T, stringsAsFactors = F)
	
	
	del = which(tmp$MarkerID < 10 | tmp$MarkerID > max(mrk$MarkerID) - 10)
	if (length(del) > 0) {
		tmp = tmp[-del, ]
	}
	del = which(tmp$LHS < 10)
	if (length(del) > 0) {
		tmp = tmp[-del, ]
	}
	del = which(tmp$RHS > max(tmp$RHS) - 10)
	if (length(del) > 0) {
		tmp = tmp[-del, ]
	}
	
	tmp$foc = mrk$Position[ tmp$MarkerID+1 ]
	if (adjust.breakpts) {
		tmp$lhs = mrk$Position[ tmp$LHS ]
		tmp$rhs = mrk$Position[ tmp$RHS+2 ]
		tmp$g.lhs = mrk$GenDist[ tmp$LHS ]
		tmp$g.rhs = mrk$GenDist[ tmp$RHS+2 ]
	} else {
		tmp$lhs = mrk$Position[ tmp$LHS+1 ]
		tmp$rhs = mrk$Position[ tmp$RHS+1 ]
		tmp$g.lhs = mrk$GenDist[ tmp$LHS+1 ]
		tmp$g.rhs = mrk$GenDist[ tmp$RHS+1 ]
	}
	
	
	data = rbind(data, tmp)
}

save(data, file = "__plotdata.RData")
# load("__plotdata.RData")


data$len = (data$rhs - data$lhs) + 1
data$g.len = (data$g.rhs - data$g.lhs)


tmp = lapply(split(data, data$type), function(x) {
	k = sprintf("%s %d %d %d %d", x$chr, x$SampleID0, x$SampleID1, x$LHS, x$RHS)
	del = which(duplicated(k))
	if (length(del)) {
		x = x[-del,]
	}
	x
})
tmp = lapply(tmp, as.data.table)
tmp = rbindlist(tmp)

table(tmp$type)

se <- function(x) sqrt(var(x)/length(x))

x = by(tmp, list(tmp$chr, tmp$type), function(x) median(x$len/1e6))
y = by(tmp, list(tmp$chr, tmp$type), function(x) median(x$g.len))
z = by(tmp, list(tmp$chr, tmp$type), function(x) mean(x$g.len / (x$len/1e6)))
zz = by(tmp, list(tmp$chr, tmp$type), function(x) se(x$g.len / (x$len/1e6)))

x = array(x, dim(x), dimnames(x))
x = x[order(as.numeric(rownames(x))),]

y = array(y, dim(y), dimnames(y))
y = y[order(as.numeric(rownames(y))),]

z = array(z, dim(z), dimnames(z))
z = z[order(as.numeric(rownames(z))),]

zz = array(zz, dim(zz), dimnames(zz))
zz = zz[order(as.numeric(rownames(zz))),]

q = cbind(x, y, z[,1], zz[,1], z[,2], zz[,2])
colnames(q) = NULL
round(q, 3)



x = which(data$chr == "20")
d = data[x, ]


table(d$type)

d$type[grep("FGT", d$type)] = "(b) FGT, phased haplotypes"
d$type[grep("DGT", d$type)] = "(c) DGT, genotypes"



if (unique.segments) {
	splt = split(d, d$type)
	
	print(range(sapply(splt, function(x) table(x$Fk))))
	
	splt = lapply(splt, function(x) {
		del = sprintf("%d %d %d %d %d", x$Fk, x$SampleID0, x$SampleID1, x$LHS, x$RHS) # per fk !!!
		del = which(duplicated(del))
		x[-del, ]
	})
	
	print(range(sapply(splt, function(x) table(x$Fk))))
	
	d = rbindlist(splt)
	
	table(d$type)
}



if (subsample.lengths) {
	table(d$Fk)
	
	tmp = lapply(split(d, list(d$type, d$Fk)), function(x) {
		if (nrow(x) < 10000) return(x)
		x[sample(1:nrow(x), 10000),]
	})
	d = rbindlist(tmp)
	table(d$Fk)
}




alt = data.table(x = as.numeric(factor(2:25))[seq(1,24, by=2)], ymin = 1e-8, ymax = 100)




p = lapply(split(d, list(d$type, d$Fk)), function(x) {
	tmp = boxplot.stats(x$len)$stats
	data.table(type = x$type[1],
						 fk = x$Fk[1],
						 lower = tmp[2] / 1e6,
						 median = tmp[3] / 1e6,
						 upper = tmp[4] / 1e6)
})
p = rbindlist(p)

gg = ggplot(p) + 
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 2, alpha = 1, position = position_dodge(width = 3/5)) +
	geom_rect(data = alt, aes(xmin = x-0.5, xmax = x+0.5, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1) +
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 2, alpha = 1, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "white", alpha = 1/1, size = 2, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, colour = type), size = 1, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "black", alpha = 0.5, size = 1, position = position_dodge(width = 3/5)) +
	scale_y_log10(expand=c(0.025, 0), breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20), 
								labels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20)), 
								minor_breaks = c(seq(0.01, 0.095, by=0.005), seq(0.1, 0.95, by=0.05), seq(1, 10, by=0.5), seq(10, 100, by=5))) +
	scale_x_discrete(expand = c(0,0)) +
	scale_colour_manual(values = c("seagreen3", "goldenrod")) +
	coord_cartesian(ylim=c(0.038, 9.2), xlim = c(0.5, 24.5)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		panel.grid = element_line(colour = "grey90", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey90", size=1/3),
		panel.grid.major.y = element_line(colour = "grey90", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		legend.justification = c(1, 1), legend.position = c(1, 1), legend.margin = unit(0.1, "cm"), legend.background = element_rect(fill = "white", colour = "black", size = 0.25), 
		#legend.position = "none",
		legend.direction = "vertical",
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Physical length (Mb)")

gg

ggsave(gg, filename = "_boxplot.1kg20.phylen.pdf", height = 4.5, width = 9)



p = lapply(split(d, list(d$type, d$Fk)), function(x) {
	tmp = boxplot.stats(x$g.len)$stats
	data.table(type = x$type[1],
						 fk = x$Fk[1],
						 lower = tmp[2],
						 median = tmp[3],
						 upper = tmp[4])
})
p = rbindlist(p)

gg = ggplot(p) + 
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 2, alpha = 1, position = position_dodge(width = 3/5)) +
	geom_rect(data = alt, aes(xmin = x-0.5, xmax = x+0.5, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1) +
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 2, alpha = 1, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "white", alpha = 1/1, size = 2, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, colour = type), size = 1, position = position_dodge(width = 3/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "black", alpha = 0.5, size = 1, position = position_dodge(width = 3/5)) +
	scale_y_log10(expand=c(0.025, 0), breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20), 
								labels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20)), 
								minor_breaks = c(seq(0.01, 0.095, by=0.005), seq(0.1, 0.95, by=0.05), seq(1, 10, by=0.5), seq(10, 100, by=5))) +
	scale_x_discrete(expand = c(0,0)) +
	scale_colour_manual(values = c("seagreen3", "goldenrod")) +
	coord_cartesian(ylim=c(0.038, 9.2), xlim = c(0.5, 24.5)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		panel.grid = element_line(colour = "grey90", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey90", size=1/3),
		panel.grid.major.y = element_line(colour = "grey90", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		#legend.justification = c(1, 1), legend.position = c(1, 1), legend.margin = unit(0.1, "cm"), legend.background = element_rect(fill = "white", colour = "black", size = 0.25), 
		legend.position = "none",
		#legend.direction = "vertical",
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Genetic length (cM)")

gg

ggsave(gg, filename = "_boxplot.1kg20.genlen.pdf", height = 4.5, width = 9)

















###
stop()

a = sort(unique(d$fk))
b = (a / 5000) * 100
frq = b
names(frq) = a

x = ((a %% 5 == 0))
frq[ !x ] = ""

gg = ggplot(p) + 
	facet_grid(tag~.) + 
	geom_rect(data = tru, aes(xmin = as.numeric(factor(fk)) - 0.45, xmax = as.numeric(factor(fk)) + 0.45, ymin = lower, ymax = upper), fill = "blue", alpha = 1/2) +
	geom_segment(data = tru, aes(x = as.numeric(factor(fk)) - 0.45, xend = as.numeric(factor(fk)) + 0.45, y = median, yend = median), colour = "blue", size = 1, alpha = 0.5) +
	geom_boxplot(aes(x = factor(fk), y = len), fill = "white", outlier.shape = NA, width = 0.75) +
	geom_segment(data = tru, aes(x = as.numeric(factor(fk)) - 0.45, xend = as.numeric(factor(fk)) + 0.45, y = median, yend = median), colour = "blue", size = 1, alpha = 0.5) +
	#geom_point(data = tru, aes(x = factor(fk), y = md), colour = "white", size = 2.5, alpha = 0.5) +
	#geom_point(data = tru, aes(x = factor(fk), y = md), colour = "blue", size = 1) +
	scale_y_continuous(expand = c(0.025, 0), breaks = seq(0, 3e7, by = 0.2e7), labels = sprintf("%.1f", seq(0, 3e7, by = 0.2e7) / 1e6), minor_breaks = seq(0, 3e7, by = 1e6)) +
	scale_x_discrete(labels = frq) +
	coord_cartesian(ylim=c(0, 1.5e7)) +
	theme_few() +
	theme(aspect.ratio = 1/2,
				panel.grid = element_line(colour = "grey80", size=0.5),
				panel.grid.minor.y = element_line(colour = "grey80", size=0.5),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size=0.5),
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (Mb)")

gg

ggsave(gg, filename = "__boxplot.naive.tru.pdf", height = 10, width = 12)

ggsave(gg, filename = "__boxplot.naive.err.pdf", height = 10, width = 12)



### stats


















