#
# generate IBD region plot
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


prefix = "history.generror_1000g"
truth.file = "../truth/ibd_mrca.truth.RData"

load(truth.file)


load(sprintf("%s.RData", prefix))

G = load.bigmatrix(sprintf("%s.G", prefix))
H = load.bigmatrix(sprintf("%s.H", prefix))
P = load.bigmatrix(sprintf("%s.P", prefix))

AF = rep(0, nrow(H)); for (i in 1:ncol(H)) { AF = AF + H[, i] }
AF = AF / ncol(H)



detect.shap.gen = function(idx, g0, g1, G) {
	g0 = G[, g0]
	g1 = G[, g1]
	
	if (g0[idx] != 1 || g1[idx] != 1) {
		return(NULL)
	}
	
	brk = abs(g0 - g1)
	
	which(brk == 2)
}


detect.shap.hap = function(idx, h0, h1, H) {
	h0m = H[, h0]
	h1m = H[, h1]
	
	h0p = if (h0 %% 2 == 0) H[, h0 - 1] else H[, h0 + 1]
	h1p = if (h1 %% 2 == 0) H[, h1 - 1] else H[, h1 + 1]
	
	if (h0m[idx] != h1m[idx] ||
			h0m[idx] == h0p[idx] ||
			h1m[idx] == h1p[idx]) {
		return(NULL)
	}
	
	brk.m = abs(h0m - h1m)
	brk.p = abs(h0p - h1p)
	
	intersect(which(brk.m == 1), which(brk.p == 1))
}




fki = FKI[sample(which(FKI$n.sharer == 8), 1), ]

g.sharer = as.numeric(strsplit(fki$g.sharer, "|", T)[[1]])
g.sharer = combn(g.sharer, 2, simplify = T)

h.sharer = as.numeric(strsplit(fki$h.sharer, "|", T)[[1]])
h.sharer = combn(h.sharer, 2, simplify = T)

p.sharer = as.numeric(strsplit(fki$p.sharer, "|", T)[[1]])
p.sharer = combn(p.sharer, 2, simplify = T)



prs = list()
foc = list()
rar = list()
tru = list()
box = list()
g.brk = list()
h.brk = list()
p.brk = list()
for (i in 1:choose(fki$n.sharer, 2)) {
	z = as.character(i)
	pair = i * 4
	
	g = g.sharer[, i]
	h = h.sharer[, i]
	p = p.sharer[, i]
	
	prs[[z]] = data.table(pair = pair - 4, indv1 = as.character(g[1]), indv2 = as.character(g[2]))
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", g[1]), FKI$g.sharer)
	a = a[sapply(strsplit(FKI$g.sharer[a], "|", T), function(x, z) { z %in% x }, g[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", g[2]), FKI$g.sharer)
	b = b[sapply(strsplit(FKI$g.sharer[b], "|", T), function(x, z) { z %in% x }, g[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		foc[[paste(z, "G")]] = data.table(pair = pair - 3, foc = FKI$pos[x])
	}
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", h[1]), FKI$h.sharer)
	a = a[sapply(strsplit(FKI$h.sharer[a], "|", T), function(x, z) { z %in% x }, h[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", h[2]), FKI$h.sharer)
	b = b[sapply(strsplit(FKI$h.sharer[b], "|", T), function(x, z) { z %in% x }, h[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		foc[[paste(z, "H")]] = data.table(pair = pair - 2, foc = FKI$pos[x])
	}
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", p[1]), FKI$p.sharer)
	a = a[sapply(strsplit(FKI$p.sharer[a], "|", T), function(x, z) { z %in% x }, p[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", p[2]), FKI$p.sharer)
	b = b[sapply(strsplit(FKI$p.sharer[b], "|", T), function(x, z) { z %in% x }, p[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		foc[[paste(z, "P")]] = data.table(pair = pair - 1, foc = FKI$pos[x])
	}
	
	x = which(truth$h0 == h[1] & truth$h1 == h[2])
	if (length(x) != 0) {
		rar[[z]] = data.table(pair = pair - 2, foc = truth$position[x])
		tru[[z]] = data.table(pair = pair - 2, tru = c(truth$brk.lhs[x], truth$brk.rhs[x]))
		box[[z]] = data.table(pair = pair - 2, lhs = truth$brk.lhs[x], rhs = truth$brk.rhs[x])
	}
	
	brk.g = detect.shap.gen(fki$index, g[1], g[2], G)
	brk.h = detect.shap.hap(fki$index, h[1], h[2], H)
	brk.p = detect.shap.hap(fki$index, p[1], p[2], P)
	
	g.brk[[z]] = data.table(pair = pair - 3, brk = POS[brk.g])
	h.brk[[z]] = data.table(pair = pair - 2, brk = POS[brk.h])
	p.brk[[z]] = data.table(pair = pair - 1, brk = POS[brk.p])
}
prs = rbindlist(prs)
foc = rbindlist(foc)
rar = rbindlist(rar)
tru = rbindlist(tru)
box = rbindlist(box)
g.brk = rbindlist(g.brk)
h.brk = rbindlist(h.brk)
p.brk = rbindlist(p.brk)


tru = unique(tru)
if (any(is.na(tru$tru))) {
	del = which(is.na(tru$tru))
	tru = tru[-del, ]
}


xlim = range(POS) / 1e6
step = (diff(xlim)) * 0.01

nil = which(is.na(box$lhs)); if (length(nil) > 0) box$lhs[nil] = min(POS)
nil = which(is.na(box$rhs)); if (length(nil) > 0) box$rhs[nil] = max(POS)

ggplot() + 
	geom_vline(xintercept = fki$pos / 1e6) +
	geom_vline(xintercept = xlim, colour = "grey", size = 1) +
	geom_linerange(data = g.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darksalmon") + 
	geom_linerange(data = h.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darkseagreen") + 
	geom_linerange(data = p.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darkslategray3") + 
	geom_linerange(data = rar, aes(x = foc / 1e6, y = pair, ymin = pair - 1.75, ymax = pair + 1.75), colour = "orange") +
	geom_linerange(data = foc, aes(x = foc / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "red") +
	geom_rect(data = box, aes(xmin = lhs / 1e6, xmax = rhs / 1e6, ymin = pair - 1.75, ymax = pair + 1.75), colour = "blue", fill = NA) +
	geom_text(data = prs, aes(x = min(xlim), y = pair, label = indv1, colour = indv1), nudge_x = step * 2, show.legend = F) +
	geom_text(data = prs, aes(x = min(xlim), y = pair, label = indv2, colour = indv2), nudge_x = step * 6, show.legend = F) +
	scale_y_reverse() +
	scale_x_continuous(expand = c(0, 0)) +
	coord_cartesian(xlim = xlim) +
	theme_classic() +
	theme(axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.line.y = element_blank(),
				axis.line.x = element_line(),
				axis.title.y = element_blank()) +
	xlab("Position (Mb)")

ggsave(filename = sprintf("_plot_sharingpairs_truth.%d.chrom.full.pdf", fki$index), width = 11.69 * 1.5, height = 8.27 * 1.5)


xlim = c(max(fki$pos - 5e6, min(POS)), min(fki$pos + 5e6, max(POS)))  / 1e6
step = (diff(xlim)) * 0.01

nil = which(is.na(box$lhs)); if (length(nil) > 0) box$lhs[nil] = min(POS)
nil = which(is.na(box$rhs)); if (length(nil) > 0) box$rhs[nil] = max(POS)

ggplot() + 
	geom_vline(xintercept = fki$pos / 1e6) +
	geom_vline(xintercept = xlim, colour = "grey", size = 1) +
	geom_linerange(data = g.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darksalmon") + 
	geom_linerange(data = h.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darkseagreen") + 
	geom_linerange(data = p.brk, aes(x = brk / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "darkslategray3") + 
	geom_linerange(data = rar, aes(x = foc / 1e6, y = pair, ymin = pair - 1.75, ymax = pair + 1.75), colour = "orange") +
	geom_linerange(data = foc, aes(x = foc / 1e6, y = pair, ymin = pair - 0.45, ymax = pair + 0.45), colour = "red") +
	geom_rect(data = box, aes(xmin = lhs / 1e6, xmax = rhs / 1e6, ymin = pair - 1.75, ymax = pair + 1.75), colour = "blue", fill = NA) +
	geom_text(data = prs, aes(x = min(xlim), y = pair, label = indv1, colour = indv1), nudge_x = step * 2, show.legend = F) +
	geom_text(data = prs, aes(x = min(xlim), y = pair, label = indv2, colour = indv2), nudge_x = step * 6, show.legend = F) +
	scale_y_reverse() +
	scale_x_continuous(expand = c(0, 0)) +
	coord_cartesian(xlim = xlim) +
	theme_classic() +
	theme(axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.line.y = element_blank(),
				axis.line.x = element_line(),
				axis.title.y = element_blank()) +
	xlab("Position (Mb)")

ggsave(filename = sprintf("_plot_sharingpairs_truth.%d.chrom.zoom.pdf", fki$index), width = 11.69 * 1.5, height = 8.27 * 1.5)



d = list()
f = list()
m = list()
for (i in 1:choose(fki$n.sharer, 2)) {
	z = as.character(i)
	pair = i
	
	g = g.sharer[, i]
	h = h.sharer[, i]
	p = p.sharer[, i]
	
	x = which(truth$position == fki$pos & truth$h0 == h[1] & truth$h1 == h[2])
	
	if (length(x) != 1 || is.na(truth$wall.lhs[x]) || is.na(truth$wall.rhs[x])) {
		next
	}
	
	x = data.table(pair = pair, 
								 indv1 = as.character(g[1]), 
								 indv2 = as.character(g[2]),
								 foc = truth$position[x],
								 lhs = truth$brk.lhs[x], 
								 rhs = truth$brk.rhs[x])
	
	brk.g = detect.shap.gen(fki$index, g[1], g[2], G)
	brk.h = detect.shap.hap(fki$index, h[1], h[2], H)
	brk.p = detect.shap.hap(fki$index, p[1], p[2], P)
	
	d[[z]] = rbind(cbind(x, brk = POS[brk.g], maf = MAF[brk.g], af = AF[brk.g], type = "Genotype"),
								 cbind(x, brk = POS[brk.h], maf = MAF[brk.h], af = AF[brk.h], type = "Haplotype"),
								 cbind(x, brk = POS[brk.p], maf = MAF[brk.p], af = AF[brk.p], type = "Phased haplotype"))
	
	
	x = which(truth$h0 == h[1] & truth$h1 == h[2])
	if (length(x) != 0) {
		m[[z]] = data.table(pair = pair, foc = truth$position[x])
	}
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", g[1]), FKI$g.sharer)
	a = a[sapply(strsplit(FKI$g.sharer[a], "|", T), function(x, z) { z %in% x }, g[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", g[2]), FKI$g.sharer)
	b = b[sapply(strsplit(FKI$g.sharer[b], "|", T), function(x, z) { z %in% x }, g[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		f[[paste(z, "G")]] = data.table(pair = pair, foc = FKI$pos[x], type = "Genotype")
	}
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", h[1]), FKI$h.sharer)
	a = a[sapply(strsplit(FKI$h.sharer[a], "|", T), function(x, z) { z %in% x }, h[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", h[2]), FKI$h.sharer)
	b = b[sapply(strsplit(FKI$h.sharer[b], "|", T), function(x, z) { z %in% x }, h[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		f[[paste(z, "H")]] = data.table(pair = pair, foc = FKI$pos[x], type = "Haplotype")
	}
	
	a = grep(sprintf("[^0-9]*%d[^0-9]*", p[1]), FKI$p.sharer)
	a = a[sapply(strsplit(FKI$p.sharer[a], "|", T), function(x, z) { z %in% x }, p[1])]
	b = grep(sprintf("[^0-9]*%d[^0-9]*", p[2]), FKI$p.sharer)
	b = b[sapply(strsplit(FKI$p.sharer[b], "|", T), function(x, z) { z %in% x }, p[2])]
	x = intersect(a, b)
	if (length(x) != 0) {
		f[[paste(z, "P")]] = data.table(pair = pair, foc = FKI$pos[x], type = "Phased haplotype")
	}
	
}
d = rbindlist(d)
f = rbindlist(f)
m = rbindlist(m)

d$foc = d$foc / 1e6
d$lhs = d$lhs / 1e6
d$rhs = d$rhs / 1e6
d$brk = d$brk / 1e6
f$foc = f$foc / 1e6
m$foc = m$foc / 1e6


d = split(d, d$pair)
f = split(f, f$pair)
m = split(m, m$pair)


for (tag in names(d)) {
	
	p = d[[tag]]
	j = f[[tag]]
	w = m[[tag]]
	
	
	plhs = ifelse(is.na(p$lhs[1]), min(p$brk), p$lhs[1])
	prhs = ifelse(is.na(p$rhs[1]), max(p$brk), p$rhs[1])
	
	dist = prhs - plhs
	step = dist * 0.25
	xlim = c(max(c(plhs - step, min(p$brk))), min(c(prhs + step, max(p$brk))))
	
	
	t = split(p, p$type)
	t = lapply(t, function(x, n) {
		ins = length(which(x$brk >  x$lhs & x$brk <  x$rhs)) / n
		out = length(which(x$brk <= x$lhs | x$brk >= x$rhs)) / n
		data.table(type = x$type[1], 
							 ins = sprintf("Breakpoints inside:  %.4f%%", ins * 100), 
							 out = sprintf("Breakpoints outside:  %.4f%%", out * 100),
							 all = sprintf("Breakpoints overall:  %.4f%%", (ins + out) * 100))
	}, length(POS))
	t = rbindlist(t)
	
	
	p = split(p, p$type)
	p = lapply(p, function(x) {
		l = which(x$brk < x$foc)
		r = which(x$brk > x$foc)
		x$dist = c(  c(x$brk[l[-1]], x$foc[1]) - x$brk[l] ,
								 (c(x$foc[1], x$brk[r[-(length(r))]]) - x$brk[r]) * -1  )
		x$dist = x$dist / max(x$dist)
		x
	})
	p = rbindlist(p)
	
	
	
	gg = ggplot(data = p) +
		facet_grid(type~.) +
		#geom_vline(data = w, aes(xintercept = foc), colour = "orange") +
		#geom_vline(data = j, aes(xintercept = foc), colour = "red", linetype="55") +
		geom_linerange(aes(x = brk, y = 0.5, ymin = 0, ymax = 1), colour = "grey60") + 
		geom_point(aes(x = brk, y = af), shape = 1) +
		geom_vline(xintercept = p$foc[1], colour = "red") +
		geom_vline(xintercept = c(plhs, prhs), colour = "blue") +
		geom_label(data = t, aes(x = plhs + (dist / 5), y = 1.05, label = ins)) +
		geom_label(data = t, aes(x = prhs - (dist / 5), y = 1.05, label = out)) +
		geom_label(data = t, aes(x = plhs + (dist / 2), y = 1.05, label = all)) +
		coord_cartesian(xlim = xlim) +
		theme_few() +
		xlab("Position (Mb)") +
		ylab("Allele frequency") +
		ggtitle(sprintf("Breakpoint allele frequency    //    Individuals:  %s  %s    //    True length:  %.1fMb", p$indv1[1], p$indv2[1], round(dist, 1)))
	
	ggsave(plot = gg, filename = sprintf("_plot_sharingpairs_truth.%d.pair-%s-%s.freq.pdf", fki$index, p$indv1[1], p$indv2[1]), width = 11.69 * 1.5, height = 8.27 * 1.5)
	
	gg = ggplot(data = p) +
		facet_grid(type~.) +
		#geom_vline(data = w, aes(xintercept = foc), colour = "orange") +
		#geom_vline(data = j, aes(xintercept = foc), colour = "red", linetype="55") +
		geom_linerange(aes(x = brk, y = 0.5, ymin = 0, ymax = 1), colour = "grey60") + 
		geom_point(aes(x = brk, y = dist), shape = 1) +
		geom_vline(xintercept = p$foc[1], colour = "red") +
		geom_vline(xintercept = c(plhs, prhs), colour = "blue") +
		geom_label(data = t, aes(x = plhs + (dist / 5), y = 1.05, label = ins)) +
		geom_label(data = t, aes(x = prhs - (dist / 5), y = 1.05, label = out)) +
		geom_label(data = t, aes(x = plhs + (dist / 2), y = 1.05, label = all)) +
		coord_cartesian(xlim = xlim) +
		theme_few() +
		xlab("Position (Mb)") +
		ylab("Breakpoint distance") +
		ggtitle(sprintf("Relative breakpoint distance    //    Individuals:  %s  %s    //    True length:  %.1fMb", p$indv1[1], p$indv2[1], round(dist, 1)))
	
	ggsave(plot = gg, filename = sprintf("_plot_sharingpairs_truth.%d.pair-%s-%s.dist.pdf", fki$index, p$indv1[1], p$indv2[1]), width = 11.69 * 1.5, height = 8.27 * 1.5)
	
	p = split(p, p$type)
	p = lapply(p, function(x) {
		l = which(x$brk < x$foc)
		r = which(x$brk > x$foc)
		x$dist = c(  c(x$brk[l[-1]], x$foc[1]) - x$brk[l] ,
								 (c(x$foc[1], x$brk[r[-(length(r))]]) - x$brk[r]) * -1  ) * 1e6
		
		ins = which(x$brk >  x$lhs & x$brk <  x$rhs)
		x$region = "Outside"
		x$region[ins] = "Inside"
		x
	})
	p = rbindlist(p)
	
	
	gg = ggplot(data = p) +
		facet_grid(type~.) +
		geom_density(aes(log10(dist+1), fill = region, colour = region), alpha = 0.5) +
		scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
		scale_y_continuous(expand = c(0, 0), breaks = seq(0, 8, by=2)/10) +
		coord_cartesian(ylim = c(0, 1)) +
		theme_base() +
		theme(legend.position = "top",
					legend.title = element_blank()) +
		xlab("Relative breakpoint distance, log10(bp)") +
		ylab("Density") +
		ggtitle(sprintf("Density of relative breakpoint distance    //    Individuals:  %s  %s    //    True length:  %.1fMb", p$indv1[1], p$indv2[1], round(dist, 1)))
	
	ggsave(plot = gg, filename = sprintf("_plot_sharingpairs_truth.%d.pair-%s-%s.distdens.pdf", fki$index, p$indv1[1], p$indv2[1]), width = 11.69 * 1.5, height = 8.27 * 1.5)
}



