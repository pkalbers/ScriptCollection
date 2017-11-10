

library(data.table)
library(ggplot2)
library(ggthemes)


unique.segments = T
equal.nsegments = F
adjust.breakpts = T
adjust.breakpts.truth = T
matched.segments = T
subsample.lengths = T



prefix = "OutOfAfricaHapMap20.GenErr_1000G" #args[1]

load(sprintf("data.%s.RData", prefix))

load(sprintf("result.beagle.%s.RData", prefix))



M = fread("~/Research/rvage/africa/truH.marker.txt", header = T, stringsAsFactors = F)

pos = M$Position
names(pos) = as.character(M$MarkerID)

gen = M$GenDist
names(gen) = as.character(M$MarkerID)




load("~/Research/DetectIBD/result.truth.local.RData")

del = which(truth$wall)
if (length(del) > 0)
	truth = truth[-del,]


if (adjust.breakpts.truth) {
	truth$length = (pos[ (truth$rhs.index + 1) ] - pos[ (truth$lhs.index - 1) ]) + 1
	truth$g.length = (gen[ (truth$rhs.index + 2) ] - gen[ (truth$lhs.index) ])
} else {
	truth$length = (pos[ (truth$rhs.index) ] - pos[ (truth$lhs.index) ]) + 1
	truth$g.length = (gen[ (truth$rhs.index + 1) ] - gen[ (truth$lhs.index + 1) ])
}
del = which(is.na(truth$length))
if (length(del) > 0)
	truth = truth[-del,]
del = which(is.na(truth$g.length))
if (length(del) > 0)
	truth = truth[-del,]


TRUTH = truth  #  truth = TRUTH







del = which(match.H$beagle.lhs.index < 10 | match.H$beagle.rhs.index > length(pos) - 10)
if (length(del) > 0) {
	match.H = match.H[-del, ]
}
del = which(match.P$beagle.lhs.index < 10 | match.P$beagle.rhs.index > length(pos) - 10)
if (length(del) > 0) {
	match.P = match.P[-del, ]
}

d = rbind(data.table(fk = match.H$fk, 
										 id0 = match.H$g0, 
										 id1 = match.H$g1, 
										 foc = match.H$index,
										 lhs = match.H$beagle.lhs.index-1,
										 rhs = match.H$beagle.rhs.index+1,
										 type = "(a) Beagle IBD, true haplotypes",   
										 length   = pos[match.H$beagle.rhs.index+1] - pos[match.H$beagle.lhs.index-1],
										 g.length = gen[match.H$beagle.rhs.index+1] - gen[match.H$beagle.lhs.index-1]
										 ),
					data.table(fk = match.P$fk, 
										 id0 = match.P$g0, 
										 id1 = match.P$g1, 
										 foc = match.P$index,
										 lhs = match.P$beagle.lhs.index-1,
										 rhs = match.P$beagle.rhs.index+1,
										 type = "(b) Beagle IBD, phased haplotypes",   
										 length   = pos[match.P$beagle.rhs.index+1] - pos[match.P$beagle.lhs.index-1],
										 g.length = gen[match.P$beagle.rhs.index+1] - gen[match.P$beagle.lhs.index-1]
					))





table(d$type)

del = which(is.na(d$rhs) | is.na(d$lhs)); if (length(del) > 0) d = d[-del,]
del = which(is.na(d$g.len) | is.na(d$g.len)); if (length(del) > 0) d = d[-del,]





if (unique.segments) {
	splt = split(d, d$type)
	
	splt = lapply(splt, function(x) {
		del = sprintf("%d %d %d %d %d", x$fk, x$id0, x$id1, x$lhs, x$rhs) # per fk !!!
		del = which(duplicated(del))
		x[-del, ]
	})
	
	d = rbindlist(splt)
	
	
	del = sprintf("%d %d %d %d %d", truth$fk, truth$g0, truth$g1, truth$lhs.index, truth$rhs.index) # per fk !!!
	del = which(duplicated(del))
	truth = truth[-del, ]
}



if (subsample.lengths) {
	tmp = lapply(split(truth, truth$fk), function(x) {
		if (nrow(x) < 10000) return(x)
		x[sample(1:nrow(x), 10000),]
	})
	truth = rbindlist(tmp)
	table(truth$fk)
}

if (subsample.lengths) {
	tmp = lapply(split(d, list(d$type, d$fk)), function(x) {
		if (nrow(x) < 10000) return(x)
		x[sample(1:nrow(x), 10000),]
	})
	d = rbindlist(tmp)
	table(d$fk)
}




alt = data.table(x = as.numeric(factor(2:25))[seq(1,24, by=2)], ymin = 1e-8, ymax = 100)




tru = lapply(split(truth, truth$fk), function(x) {
	tmp = boxplot.stats(x$length)$stats
	data.table(fk = x$fk[1],
						 lower = tmp[2] / 1e6,
						 median = tmp[3] / 1e6,
						 upper = tmp[4] / 1e6)
})
tru = rbindlist(tru)



det = lapply(split(d, list(d$type, d$fk)), function(x) {
	tmp = boxplot.stats(x$length)$stats
	data.table(type = x$type[1],
						 fk = x$fk[1],
						 lower = tmp[2] / 1e6,
						 median = tmp[3] / 1e6,
						 upper = tmp[4] / 1e6)
})
det = rbindlist(det)

tru = cbind(type = " True IBD", tru)
p = rbind(det, tru)





gg = ggplot(p) + 
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 1.75, alpha = 1, position = position_dodge(width = 4/5)) +
	geom_rect(data = alt, aes(xmin = x-0.5, xmax = x+0.5, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1) +
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 1.75, alpha = 1, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "white", alpha = 1/1, size = 1.5, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, colour = type), size = 1/2, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "black", alpha = 0.25, size = 1/2, position = position_dodge(width = 4/5)) +
	scale_y_log10(expand=c(0.025, 0), breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20), 
								labels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20)), 
								minor_breaks = c(seq(0.01, 0.095, by=0.005), seq(0.1, 0.95, by=0.05), seq(1, 10, by=0.5), seq(10, 100, by=5))) +
	scale_x_discrete(expand = c(0,0)) +
	scale_colour_manual(values = c("black", "purple", "seagreen3", "goldenrod")) + #"plum3", "goldenrod")) +
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
		#legend.position = "bottom",
		legend.direction = "vertical",
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Physical length (Mb)")

gg

ggsave(gg, filename = "_boxplot.phylen.beagle.tru.pdf", height = 4.5, width = 9)
ggsave(gg, filename = "_boxplot.phylen.beagle.err.pdf", height = 4.5, width = 9)




tru = lapply(split(truth, truth$fk), function(x) {
	tmp = boxplot.stats(x$g.length)$stats
	data.table(fk = x$fk[1],
						 lower = tmp[2],
						 median = tmp[3],
						 upper = tmp[4])
})
tru = rbindlist(tru)



det = lapply(split(d, list(d$type, d$fk)), function(x) {
	tmp = boxplot.stats(x$g.length)$stats
	data.table(type = x$type[1],
						 fk = x$fk[1],
						 lower = tmp[2],
						 median = tmp[3],
						 upper = tmp[4])
})
det = rbindlist(det)

tru = cbind(type = " True IBD", tru)
p = rbind(det, tru)


gg = ggplot(p) + 
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 1.75, alpha = 1, position = position_dodge(width = 4/5)) +
	geom_rect(data = alt, aes(xmin = x-0.5, xmax = x+0.5, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1) +
	geom_linerange(aes(x = (factor(fk)), ymin = lower, ymax = upper, colour = type), size = 1.75, alpha = 1, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "white", alpha = 1/1, size = 1.5, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, colour = type), size = 1/2, position = position_dodge(width = 4/5)) +
	geom_point(aes(x = (factor(fk)), y = median, group = type), colour = "black", alpha = 0.25, size = 1/2, position = position_dodge(width = 4/5)) +
	scale_y_log10(expand=c(0.025, 0), breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20), 
								labels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 20)), 
								minor_breaks = c(seq(0.01, 0.095, by=0.005), seq(0.1, 0.95, by=0.05), seq(1, 10, by=0.5), seq(10, 100, by=5))) +
	scale_x_discrete(expand = c(0,0)) +
	scale_colour_manual(values = c("black", "purple", "seagreen3", "goldenrod")) + #"plum3", "goldenrod")) +
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

ggsave(gg, filename = "_boxplot.genlen.beagle.tru.pdf", height = 4.5, width = 9)
ggsave(gg, filename = "_boxplot.genlen.beagle.err.pdf", height = 4.5, width = 9)
















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


















