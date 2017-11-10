#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)

library(rPython)

#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("data.%s.RData", prefix))
load(sprintf("result.match.%s.RData", prefix))

load(sprintf("data.%s.H.RData", prefix))
load(sprintf("data.%s.G.RData", prefix))


hist.file = "./OutOfAfricaHapMap20.hdf5"           ###### path to simulation history file
load.file = "./_exec_msprime/msprime.fetch_ibd.py" ###### path to fetch file


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





del = which(match$g.wall | match$h.wall | match$p.wall | match$true.wall)
if (length(del) > 0) {
	match = match[-del, ]
}


z = rbind(data.table(index = match$index, fk = match$fk, g0 = match$g0, g1 = match$g1,  det = match$pos - match$g.lhs, tru = match$pos - POS[match$true.lhs.idx]),
					data.table(index = match$index, fk = match$fk, g0 = match$g0, g1 = match$g1,  det = match$g.rhs - match$pos, tru = POS[match$true.rhs.idx] - match$pos))

z$diff = z$det - z$tru


### random

i = sample(which(z$fk == 4 & z$index > 500000 & z$index < 510000), 1)
fki = FKI[which(FKI$index == z$index[i]), ]


d = NULL
f = NULL
r = NULL

for (i in 1:nrow(fki)) {
	h0 = fki$h0[i]
	h1 = fki$h1[i]
	
	pair = sprintf("%d + %d", h0, h1)
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
	scale_y_reverse(breaks = cmb, labels = names(cmb)) +
	coord_cartesian(expand = F) +
	theme_few() + 
	theme(aspect.ratio = 1/n,
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

gg

ggsave(gg, filename = sprintf("_plot.full_ibd.%d.%s.map.png", fki$index[1], prefix), width = 12, height = (n + 1) * 1.25)


gg = ggplot(r) +
	geom_line(aes(x = pos, y = den, colour = pair)) +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6), expand = c(0,0)) +
	scale_color_solarized() +
	theme_few() + 
	theme(aspect.ratio = 1/n,
				legend.title = element_blank(),
				legend.position = "bottom") +
	xlab("Physical position (Mb)") +
	ylab("Normalised proportion")

gg

ggsave(gg, filename = sprintf("_plot.full_ibd.%d.%s.frq.png", fki$index[1], prefix), width = 12, height = (n + 1) * 1.25)



### under/overestimated


udr = which(z$diff < 0)
udr = z[udr, ]
udr = udr[order(udr$diff, decreasing = F), ]

ovr = which(z$diff > 0)
ovr = z[ovr, ]
ovr = ovr[order(ovr$diff, decreasing = T), ]

mid = which(z$tru > mean(z$tru) * 8 & z$tru < mean(z$tru) * 10)
mid = z[mid, ]
mid = mid[order(mid$diff, decreasing = T), ]



# testing....

i = sample(which(udr$fk == 3), 1)
i = 10
fki = FKI[which(FKI$index == udr$index[i] & FKI$g0 == udr$g0[i] & FKI$g1 == udr$g1[i]), ]
fki

i = sample(which(ovr$fk == 5), 1)
fki = FKI[which(FKI$index == ovr$index[i] & FKI$g0 == ovr$g0[i] & FKI$g1 == ovr$g1[i]), ]
fki

for (qqq in 1:20) {
i = sample(which(mid$fk == 3 & mid$index > 1.25e5 & mid$index < 2.5e5), 1)
fki = FKI[which(FKI$index == mid$index[i] & FKI$g0 == mid$g0[i] & FKI$g1 == mid$g1[i]), ]
fki


d = NULL
r = NULL
f = NULL

fgt = rep(F, nrow(H))
dgt = rep(F, nrow(G))

for (i in 1:nrow(fki)) {
	g0 = fki$g0[i]
	g1 = fki$g1[i]
	
	h0 = fki$h0[i]
	h1 = fki$h1[i]
	
	A0 = if (g0 %% 2 != 0) g0 * 2 else g0 * 2 - 1
	A1 = if (g0 %% 2 != 0) g0 * 2 - 1 else g0 * 2
	B0 = if (g1 %% 2 != 0) g1 * 2 else g1 * 2 - 1
	B1 = if (g1 %% 2 != 0) g1 * 2 - 1 else g1 * 2
	
	pair = sprintf("%d + %d", g0, g1)
	cat(pair, "\n")
	
	A0.A1 = fetch.ibd(A0, A1)
	A0.B0 = fetch.ibd(A0, B0)
	A0.B1 = fetch.ibd(A0, B1)
	A1.B0 = fetch.ibd(A1, B0)
	A1.B1 = fetch.ibd(A1, B1)
	B0.B1 = fetch.ibd(B0, B1)
	
	d = rbind(d, 
						cbind(pair = pair, cmb = "A1, A2", A0.A1),
						cbind(pair = pair, cmb = "A1, B1", A0.B0),
						cbind(pair = pair, cmb = "A1, B2", A0.B1),
						cbind(pair = pair, cmb = "A2, B1", A1.B0),
						cbind(pair = pair, cmb = "A2, B2", A1.B1),
						cbind(pair = pair, cmb = "B1, B2", B0.B1))
	
	x0 = which(FKI$h0 == A0 & FKI$h1 == A1)
	x1 = which(FKI$h0 == A0 & FKI$h1 == B0)
	x2 = which(FKI$h0 == A0 & FKI$h1 == B1)
	x3 = which(FKI$h0 == A1 & FKI$h1 == B0)
	x4 = which(FKI$h0 == A1 & FKI$h1 == B1)
	x5 = which(FKI$h0 == B0 & FKI$h1 == B1)
	
	if (length(x0) > 0) r = rbind(r, data.table(pair = pair, cmb = "A1, A2", pos = FKI$pos[x0], fk = AAC[FKI$index[x0]]))
	if (length(x1) > 0) r = rbind(r, data.table(pair = pair, cmb = "A1, B1", pos = FKI$pos[x1], fk = AAC[FKI$index[x1]]))
	if (length(x2) > 0) r = rbind(r, data.table(pair = pair, cmb = "A1, B2", pos = FKI$pos[x2], fk = AAC[FKI$index[x2]]))
	if (length(x3) > 0) r = rbind(r, data.table(pair = pair, cmb = "A2, B1", pos = FKI$pos[x3], fk = AAC[FKI$index[x3]]))
	if (length(x4) > 0) r = rbind(r, data.table(pair = pair, cmb = "A2, B2", pos = FKI$pos[x4], fk = AAC[FKI$index[x4]]))
	if (length(x5) > 0) r = rbind(r, data.table(pair = pair, cmb = "B1, B2", pos = FKI$pos[x5], fk = AAC[FKI$index[x5]]))
	
	if (A0 == h0 && A1 == h1) f = rbind(f, data.table(pair = pair, cmb = "A1, A2", pos = fki$pos[i]))
	if (A0 == h0 && B0 == h1) f = rbind(f, data.table(pair = pair, cmb = "A1, B1", pos = fki$pos[i]))
	if (A0 == h0 && B1 == h1) f = rbind(f, data.table(pair = pair, cmb = "A1, B2", pos = fki$pos[i]))
	if (A1 == h0 && B0 == h1) f = rbind(f, data.table(pair = pair, cmb = "A2, B1", pos = fki$pos[i]))
	if (A1 == h0 && B1 == h1) f = rbind(f, data.table(pair = pair, cmb = "A2, B2", pos = fki$pos[i]))
	if (B0 == h0 && B1 == h1) f = rbind(f, data.table(pair = pair, cmb = "B1, B2", pos = fki$pos[i]))
	
	
	x = matrix(NA, nrow = nrow(H), ncol = 4)
	x[, 1] = H[, A0]
	x[, 2] = H[, A1]
	x[, 3] = H[, B0]
	x[, 4] = H[, B1]
	y = apply(x, 1, sum)
	idx = fki$index[1]
	foc = x[idx, ]
	for (i in 1:nrow(H)) {
		if (y[i] != 2) next
		if (length(unique(paste(foc, x[i, ]))) == 4) fgt[i] = T
	}
	
	x = abs(G[, g0] - G[, g1])
	dgt[which(x == 2)] = T
}

x = which(d$time < 1);  if (length(x) > 0) d$time[x] = 1

tmp = unique(d$cmb)
cmb = 1:length(tmp)
names(cmb) = tmp

d$cmb = cmb[d$cmb]
r$cmb = cmb[r$cmb]
f$cmb = cmb[f$cmb]

gg = ggplot(d) + 
	geom_rect(aes(xmin = lhs.position, xmax = rhs.position, ymin = cmb - 0.5, ymax = cmb + 0.5, fill = log10(time))) + 
	geom_hline(yintercept = (0:6)+0.5, color = "white", size = 2) +
	#geom_point(data = r, aes(x = pos, y = cmb), shape = 1, colour = "white", size = 0.5) +
	geom_point(data = f, aes(x = pos, y = cmb), shape = 3, colour = "black", size = 3) +
	geom_segment(data = data.frame(x = rep(range(c(d$lhs.position, d$rhs.position)), 2), y = c(0.45, 0.45), yend = c(7, 7) - 0.45), aes(x = x, xend = x, y = y, yend = yend), size = 1.125) +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	scale_y_reverse(breaks = cmb, labels = names(cmb)) +
	scale_fill_gradientn(colors = c("orange", "darkorange", "royalblue1", "royalblue3", "black"), na.value = "black", limits=c(0, 5), breaks = (0:5), labels = sprintf("%d", 10^(0:5))) +
	coord_cartesian(expand = F) +
	theme_few() + 
	theme(#aspect.ratio = 1/3.5,
				legend.position = "top",
				legend.key.width = unit(1, "inches"),
				legend.key.height = unit(0.2, "inches"),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.title.x = element_blank(),
				axis.ticks.y = element_blank()) +
	#xlab("Physical position (Mb)") +
	ylab("") +
	labs(fill='TMRCA (generations)')

gg

ggsave(gg, filename = sprintf("_plot.fullIBD.%d.each.%s.A.pdf", fki$index[1], prefix), width = 10, height = 3.5)
ggsave(gg, filename = sprintf("_plot.fullIBD.%d.each.%s.A.png", fki$index[1], prefix), width = 10, height = 3.5)



# FGT + DGT

p = rbind(data.table(test = "FGT", pos = POS[which(fgt)]), 
					data.table(test = "DGT", pos = POS[which(dgt)]))
p$test = factor(p$test, levels = c("FGT", "DGT"))

gg = ggplot(p) + 
	facet_wrap(~test, ncol = 1, switch = 'y') +
	geom_vline(aes(xintercept = pos), size = 0.25, color = "grey50") +
	scale_x_continuous(breaks = seq(5e6, 100e6, by = 5e6), labels = sprintf("%d", seq(5e6, 100e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F) +
	theme_few() + 
	theme(#aspect.ratio = 1/15,
				strip.switch.pad.wrap = unit(0.15, "cm"),
				strip.text = element_text(face = "bold", size = 12, angle = 45),
				panel.border = element_rect(fill = NA, colour = "black"),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				#axis.text.y = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") +
	ylab("")

gg

ggsave(gg, filename = sprintf("_plot.fullIBD.%d.each.%s.B.pdf", fki$index[1], prefix), width = 10, height = 1.5)
ggsave(gg, filename = sprintf("_plot.fullIBD.%d.each.%s.B.png", fki$index[1], prefix), width = 10, height = 1.5)


}










