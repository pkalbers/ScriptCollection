#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)

library(rPython)

#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("~/Research/DetectIBD/data.%s.RData", prefix))
load(sprintf("~/Research/DetectIBD/result.match.%s.RData", prefix))

load(sprintf("~/Research/DetectIBD/data.%s.H.RData", prefix))


hist.file = "/Users/pkalbers/Research/DetectIBD/OutOfAfricaHapMap20.hdf5"           ###### path to simulation history file
load.file = "/Users/pkalbers/Research/DetectIBD/_exec_msprime/msprime.fetch_ibd.py" ###### path to fetch file


fetch.ibd = function(h0, h1, return.raw = F, history = hist.file, load = load.file) {
	python.load(load)
	ibd = python.call("MRCA", h0, h1, history)
	ibd = data.table(mrca = ibd$mrca, 
									 time = ibd$time, 
									 lhs.position = ibd$lhs.position, 
									 rhs.position = ibd$rhs.position, 
									 lhs.index = ibd$lhs.index, 
									 rhs.index = ibd$rhs.index)
	
	if (return.raw) {
		return (ibd)
	}
	
	above = which(ibd$mrca[-(nrow(ibd))] != ibd$mrca[-1])
	below = c(1, above + 1)
	above = c(above, nrow(ibd))
	
	data.table(mrca = ibd$mrca[below], 
						 time = ibd$time[below], 
						 lhs.position = ibd$lhs.position[below], 
						 rhs.position = ibd$rhs.position[above], 
						 lhs.index = ibd$lhs.index[below], 
						 rhs.index = ibd$rhs.index[above])
}



fki = FKI[which(FKI$index == (fid + 1)), ]   ############# see below for 'fid'


d = NULL
f = NULL
r = NULL

for (i in 1:nrow(fki)) {
	h0 = fki$h0[i]
	h1 = fki$h1[i]
	
	g0 = fki$g0[i]
	g1 = fki$g1[i]
	
	pair = sprintf("%d + %d", g0, g1)
	cat(pair, "\n")
	
	ibd = fetch.ibd(h0, h1)
	
	d = rbind(d, cbind(pair = pair, cmb = i, ibd))
	
	f = rbind(f, data.table(pair = pair, cmb = i, pos = fki$pos[i]))
	
	a0 = H[, h0]
	a1 = H[, h1]
	sh = which(a0 == 1 & a1 == 1)
	sh = sh[which(AAF[sh] <= 0.01)]
	sh = density(POS[sh], bw = 1e5, n = 500, from = min(POS), to = max(POS))
	
	ps = sh$x
	ds = sh$y
	ds = ds / sum(ds)
	
	r = rbind(r, data.table(pair = pair, cmb = i, pos = ps, den = ds))
}

tmp = unique(d$pair)
cmb = 1:length(tmp)
names(cmb) = tmp

n = nrow(fki)


gg = ggplot(d) + 
	geom_rect(aes(xmin = lhs.position, xmax = rhs.position, ymin = cmb - 0.5, ymax = cmb + 0.5, fill = log10(time))) + 
	geom_hline(yintercept = (0:n)+(0.5), colour = "white", size = 2) +
	geom_point(data = f, aes(x = pos, y = cmb), colour = "white", size = 1) +
	geom_segment(data = data.frame(x = rep(range(c(d$lhs.position, d$rhs.position)), 2), y = c(0.25, 0.25), yend = n + c(0.75, 0.75)), aes(x = x, xend = x, y = y, yend = yend), size = 1) +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_fill_gradient2(high = "black", mid = "royalblue3", low = "orange", na.value = "black", midpoint = 3, limits=c(0, 5), breaks = (0:5), labels = sprintf("%d", 10^(0:5))) +
	scale_y_continuous(breaks = cmb, labels = names(cmb)) +
	coord_cartesian(expand = F) +
	theme_few() + 
	theme(aspect.ratio = 1/8,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				panel.border = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				#axis.text.y = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(gg, filename = sprintf("_plot_ibd.example.%d.true.png", fki$index[1]), width = 12, height = 6)






markers = fread("truH.marker.txt", header = T, stringsAsFactors = F)

pos = markers$Position
names(pos) = as.character(markers$MarkerID)


ff = function(file, pos) {
	x = fread(file, header = T, stringsAsFactors = F)
	x$pair = sprintf("%d + %d", x$SampleID0, x$SampleID1)
	x$type = sub("^ibd_plot_(.+)\\.ibd\\.distr\\.(.+)\\.txt$", "\\1", file)
	x$mode = sub("^ibd_plot_(.+)\\.ibd\\.distr\\.(.+)\\.txt$", "\\2", file)
	x$pos_FocalID = pos[as.character(x$FocalID)]
	if (x$mode[1] == "HMM") {
		x$pos_MarkerID = pos[as.character(x$MarkerID)]
	} else {
		x$pos_BreakID = pos[as.character(x$BreakID)]
	}
	x
}

f = dir(pattern = "^ibd_plot_.+\\.ibd\\.distr\\..+\\.txt$")
p = lapply(f, ff, pos)
names(p) = sub("^ibd_plot_(.+)\\.ibd\\.distr\\.(.+)\\.txt$", "\\1.\\2", f)

fid = lapply(p, function(x) unique(x$FocalID))
fid = Reduce(intersect, fid)
fid = fid[3]

p = lapply(p, function(x, i) x[which(x$FocalID == i), ], fid)



tru.DHG = ggplot(p$tru.DHG) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_BreakID), colour = "grey50") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/25,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(tru.DHG, filename = sprintf("_plot_ibd.example.%d.tru.DHG.png", fki$index[1]), width = 12, height = 6)


err.DHG = ggplot(p$err.DHG) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_BreakID), colour = "grey50") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/25,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(err.DHG, filename = sprintf("_plot_ibd.example.%d.err.DHG.png", fki$index[1]), width = 12, height = 6)



tru.FGT = ggplot(p$tru.FGT) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_BreakID), colour = "grey50") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/25,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(tru.FGT, filename = sprintf("_plot_ibd.example.%d.tru.FGT.png", fki$index[1]), width = 12, height = 6)


err.FGT = ggplot(p$err.FGT) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_BreakID), colour = "grey50") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/25,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(err.FGT, filename = sprintf("_plot_ibd.example.%d.err.FGT.png", fki$index[1]), width = 12, height = 6)






tru.HMM.vp = ggplot(p$tru.HMM) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_MarkerID, colour = ViterbiPath)) +
	geom_line(aes(x = pos_MarkerID, y = PostProb_IBD), colour = "white", size = 1.5, alpha = 0.5) +
	geom_line(aes(x = pos_MarkerID, y = PostProb_IBD), colour = "black", size = 0.5) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_manual(values = c(IBD="orange", NON="royalblue3")) +
	coord_cartesian(expand = F, ylim = c(-0.25, 1.25), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/16,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(tru.HMM.vp, filename = sprintf("_plot_ibd.example.%d.tru.HMM.vp.png", fki$index[1]), width = 12, height = 6)


tru.HMM.pp = ggplot(p$tru.HMM) + 
	facet_grid(pair~.) +
	geom_line(aes(x = pos_MarkerID, y = PostProb_IBD), colour = "black", size = 0.5) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_manual(values = c(IBD="orange", NON="royalblue3")) +
	coord_cartesian(expand = F, ylim = c(-0.15, 1.15), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/28,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(tru.HMM.pp, filename = sprintf("_plot_ibd.example.%d.tru.HMM.pp.png", fki$index[1]), width = 12, height = 6)


rm = p$tru.HMM
rm = split(rm, rm$pair)
rm = lapply(rm, function(x) {
	x = x[order(x$pos_MarkerID), ]
	data.table(pos = rollmean(x$pos_MarkerID / 1e6, 5000) * 1e6, val = rollmean((x$ObsProb_NON), 5000) - rollmean((x$ObsProb_IBD), 5000), pair = x$pair[1])
})
rm = rbindlist(rm)

tru.HMM.obs = ggplot(rm) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos, colour = val)) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_gradient2(high = "royalblue2", low = "orange", mid = "royalblue4", na.value = "black", midpoint = mean(range(rm$val))) +
	coord_cartesian(expand = F, ylim = c(-0.15, 1.15), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/32,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(tru.HMM.obs, filename = sprintf("_plot_ibd.example.%d.tru.HMM.obs.png", fki$index[1]), width = 12, height = 6)




err.HMM.vp = ggplot(p$err.HMM) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_MarkerID, colour = ViterbiPath)) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_manual(values = c(IBD="orange", NON="royalblue3")) +
	coord_cartesian(expand = F, ylim = c(-0.15, 1.15), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/42,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(err.HMM.vp, filename = sprintf("_plot_ibd.example.%d.err.HMM.vp.png", fki$index[1]), width = 12, height = 6)


err.HMM.pp = ggplot(p$err.HMM) + 
	facet_grid(pair~.) +
	geom_line(aes(x = pos_MarkerID, y = PostProb_IBD), colour = "black", size = 0.5) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_manual(values = c(IBD="orange", NON="royalblue3")) +
	coord_cartesian(expand = F, ylim = c(-0.15, 1.15), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/28,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(err.HMM.pp, filename = sprintf("_plot_ibd.example.%d.err.HMM.pp.png", fki$index[1]), width = 12, height = 6)


rm = p$err.HMM
rm = split(rm, rm$pair)
rm = lapply(rm, function(x) {
	x = x[order(x$pos_MarkerID), ]
	data.table(pos = rollmean(x$pos_MarkerID / 1e6, 5000) * 1e6, val = rollmean((x$ObsProb_NON), 5000) - rollmean((x$ObsProb_IBD), 5000), pair = x$pair[1])
})
rm = rbindlist(rm)

err.HMM.obs = ggplot(rm) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos, colour = val)) +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=p$tru.DHG$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_colour_gradient2(high = "royalblue2", low = "orange", mid = "royalblue4", na.value = "black", midpoint = mean(range(rm$val))) +
	coord_cartesian(expand = F, ylim = c(-0.15, 1.15), xlim = range(pos)) +
	theme_few() + 
	theme(aspect.ratio = 1/42,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.ticks.x = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

ggsave(err.HMM.obs, filename = sprintf("_plot_ibd.example.%d.err.HMM.obs.png", fki$index[1]), width = 12, height = 6)














