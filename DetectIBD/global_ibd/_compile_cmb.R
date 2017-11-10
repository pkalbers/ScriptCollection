#
# time & freq dependend genotype pair emissions
#

library(data.table)

args = commandArgs(T)
prefix = args[1]

cat("Loading...")

load(sprintf("./data.%s.RData", prefix))
#load(sprintf("./data.%s.G.RData", prefix)) # genotypes
#load(sprintf("./data.%s.H.RData", prefix)) # haplotypes
load(sprintf("./data.%s.P.RData", prefix)) # phased haplotypes

cat("OK\n")



prune.ibd = function(d) {
	d$mrca = NULL
	d$lhs.position = NULL
	d$rhs.position = NULL
	d$lhs.index = d$lhs.index + 2
	d$rhs.index = d$rhs.index - 2
	l = d$rhs.index - d$lhs.index + 1
	x = which(l < 10)
	if (length(x) > 0) {
		d = d[-x,]
	}
	if (nrow(d) == 0) {
		return(NULL)
	}
	d
}


cut.time = function(t, brk = c(0, 10^(seq(0, 6, length.out = 100))), lab = sprintf("%0.4f", 10^(seq(0, 6, length.out = 100)) )) {
	as.character(cut(t, breaks = brk, labels = lab, include.lowest = T))
}

cut.freq = function(t, brk = c(0, seq(0, 1, length.out = 101)[-1]), lab = sprintf("%0.2f", seq(0, 1, length.out = 101)[-1] )) {
	as.character(cut(t, breaks = brk, labels = lab, include.lowest = T))
}


#combine.rng = function(mat, len = nrow(G)) {
#combine.rng = function(mat, len = nrow(H)) {
combine.rng = function(mat, len = nrow(P)) {
	out = rep(F, len)
	rng = apply(mat, 1, function(x) {
		x[2] : x[3]
	})
	for (i in 1:length(rng)) {
		r = rng[[i]]
		out[r] = T
	}
	out
}


make.gen.array = function() {
	tt = 10^(seq(0, 6, length.out = 100))
	ff = seq(0, 1, length.out = 101)[-1]
	out = list()
	for (t in tt) {
		tmp = list()
		for (f in ff) {
			m = rep(0, 6)
			names(m) = c("00", "01", "02", "11", "12", "22")
			tmp[[sprintf("%0.2f", f)]] = m
		}
		out[[sprintf("%0.4f", t)]] = tmp
	}
	out
}

make.hap.array = function() {
	tt = 10^(seq(0, 6, length.out = 100))
	ff = seq(0, 1, length.out = 101)[-1]
	out = list()
	for (t in tt) {
		tmp = list()
		for (f in ff) {
			m = rep(0, 3)
			names(m) = c("00", "01", "11")
			tmp[[sprintf("%0.2f", f)]] = m
		}
		out[[sprintf("%0.4f", t)]] = tmp
	}
	out
}



FRQ = cut.freq(AAF)


ibd.files = sample(dir(pattern = "^global_ibd\\..+\\.RData", path = "../", full.names = T))

for (ibd.file in ibd.files) {
	
	Sys.sleep(runif(1, min = 1, max = 5))
	
	cat(basename(ibd.file), "\n")
	
	save.file = sprintf("result_cmb_P.%s", basename(ibd.file))
	flag.file = sprintf("_flag_cmb_P.%s.txt", basename(ibd.file))
	
	if (file.exists(save.file)) next
	if (file.exists(flag.file)) next
	
	cat(".", file = flag.file)
	
	
	A0A1 = make.hap.array()
	A0B0 = make.hap.array()
	A0B1 = make.hap.array()
	A1B0 = make.hap.array()
	A1B1 = make.hap.array()
	B0B1 = make.hap.array()
	
	
	load(ibd.file)
	
	i = 0
	for (pair in names(ibd)) {
		i = i + 1
		cat(sprintf("%d: %s\n", i, pair))
		
		cmb = ibd[[pair]]
		
		cmb = lapply(cmb, prune.ibd)
		
		a = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", pair))
		b = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", pair))
		
		a0 = a * 2 - 1;  a1 = a * 2
		b0 = b * 2 - 1;  b1 = b * 2
		
		if (sample(c(0,1), 1) == 1) {
			a1 = a * 2 - 1;  a0 = a * 2
		}
		if (sample(c(0,1), 1) == 1) {
			b1 = b * 2 - 1;  b0 = b * 2
		}
		
		
		a0a1 = sprintf("%d%d", P[, a0], P[, a1]);  x = which(a0a1 == "10");  if (length(x) > 0) a0a1[x] = "01"
		a0b0 = sprintf("%d%d", P[, a0], P[, b0]);  x = which(a0b0 == "10");  if (length(x) > 0) a0b0[x] = "01"
		a0b1 = sprintf("%d%d", P[, a0], P[, b1]);  x = which(a0b1 == "10");  if (length(x) > 0) a0b1[x] = "01"
		a1b0 = sprintf("%d%d", P[, a1], P[, b0]);  x = which(a1b0 == "10");  if (length(x) > 0) a1b0[x] = "01"
		a1b1 = sprintf("%d%d", P[, a1], P[, b1]);  x = which(a1b1 == "10");  if (length(x) > 0) a1b1[x] = "01"
		b0b1 = sprintf("%d%d", P[, b0], P[, b1]);  x = which(b0b1 == "10");  if (length(x) > 0) b0b1[x] = "01"
		
		
		mat = cmb[[ sprintf("%d %d", a1, b0) ]] ### main
		
		crd = cut.time(mat$time)
		crd = split(mat, crd)
		crd = lapply(crd, combine.rng)
		crd = lapply(crd, function(x, f) {
			split((1:nrow(P))[x], f[x])
		}, FRQ)
		
		
		for (t in names(crd)) {
			for (f in names(crd[[t]])) {
				rng = crd[[t]][[f]]
				
				if (length(rng) == 0) next
				
				n = table(a0a1[rng]);  for (tag in names(n))  A0A1[[t]][[f]][tag] = A0A1[[t]][[f]][tag] + n[tag]
				n = table(a0b0[rng]);  for (tag in names(n))  A0B0[[t]][[f]][tag] = A0B0[[t]][[f]][tag] + n[tag]
				n = table(a0b1[rng]);  for (tag in names(n))  A0B1[[t]][[f]][tag] = A0B1[[t]][[f]][tag] + n[tag]
				n = table(a1b0[rng]);  for (tag in names(n))  A1B0[[t]][[f]][tag] = A1B0[[t]][[f]][tag] + n[tag]
				n = table(a1b1[rng]);  for (tag in names(n))  A1B1[[t]][[f]][tag] = A1B1[[t]][[f]][tag] + n[tag]
				n = table(b0b1[rng]);  for (tag in names(n))  B0B1[[t]][[f]][tag] = B0B1[[t]][[f]][tag] + n[tag]
			}
		}
		
	}
	
	save(A0A1, A0B0, A0B1, A1B0, A1B1, B0B1, file = save.file)
	unlink(flag.file)
}





stop()


path = "err"

a0a1 = make.hap.array()
a0b0 = make.hap.array()
a0b1 = make.hap.array()
a1b0 = make.hap.array()
a1b1 = make.hap.array()
b0b1 = make.hap.array()

res.files = dir(pattern = "^result_cmb_P\\..+", path = path, full.names = T)

for (res.file in res.files) {
	print(res.file)
	load(res.file)
	
	for (t in names(A0A1)) { for (f in names(A0A1[[t]])) { a0a1[[t]][[f]] = a0a1[[t]][[f]] + A0A1[[t]][[f]] } }
	for (t in names(A0B0)) { for (f in names(A0B0[[t]])) { a0b0[[t]][[f]] = a0b0[[t]][[f]] + A0B0[[t]][[f]] } }
	for (t in names(A0B1)) { for (f in names(A0B1[[t]])) { a0b1[[t]][[f]] = a0b1[[t]][[f]] + A0B1[[t]][[f]] } }
	for (t in names(A1B0)) { for (f in names(A1B0[[t]])) { a1b0[[t]][[f]] = a1b0[[t]][[f]] + A1B0[[t]][[f]] } }
	for (t in names(A1B1)) { for (f in names(A1B1[[t]])) { a1b1[[t]][[f]] = a1b1[[t]][[f]] + A1B1[[t]][[f]] } }
	for (t in names(B0B1)) { for (f in names(B0B1[[t]])) { b0b1[[t]][[f]] = b0b1[[t]][[f]] + B0B1[[t]][[f]] } }
	
}


compile.hap = function(tag) {
	hh = get(tag, envir = globalenv())
	for (t in names(hh)) {
		for (f in names(hh[[t]])) {
			if (sum(hh[[t]][[f]]) == 0) {
				hh[[t]][[f]] = NULL
				next
			}
			hh[[t]][[f]] = hh[[t]][[f]] / sum(hh[[t]][[f]])
			hh[[t]][[f]] = as.data.table(as.list(hh[[t]][[f]]))
			hh[[t]][[f]] = cbind(time = as.numeric(t), freq = as.numeric(f), hh[[t]][[f]])
		}
		if (length(hh[[t]]) == 0) {
			hh[[t]] = NULL
			next
		}
		hh[[t]] = rbindlist(hh[[t]])
	}
	hh = rbindlist(hh)
	hh = rbind(data.table(mode = tag, time = hh$time, freq = hh$freq, type = "00", prob = hh$`00`),
						 data.table(mode = tag, time = hh$time, freq = hh$freq, type = "01", prob = hh$`01`),
						 data.table(mode = tag, time = hh$time, freq = hh$freq, type = "11", prob = hh$`11`))
	hh[order(hh$time, decreasing = T),]
}


d = rbind(compile.hap("a0a1"),
					compile.hap("a0b0"),
					compile.hap("a0b1"),
					compile.hap("a1b0"),
					compile.hap("a1b1"),
					compile.hap("b0b1"))

del = which(d$time > 100000)
d = d[-del,]

p = ggplot(d) +
	facet_grid(mode~type) +
	geom_point(aes(x = freq - (1/200), y = prob, color = time), size = 2/3, alpha = 2/3, shape = 16) +
	scale_color_gradientn(colours = c("yellow", "orange", "red", "purple", "blue", "black"), trans = "log10", name = "TMRCA",
												breaks = 10^(0:6), labels = format(10^(0:6), trim = T, scientific = F, big.mark = ',')) +
	scale_linetype_manual(values = c("11", "22")) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				panel.background = element_rect(fill = "grey95"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("Probability")
p

ggsave(p, filename = sprintf("_plot.H.pairs.%s.pdf", path), width = 7.5, height = 10)
ggsave(p, filename = sprintf("_plot.P.pairs.%s.pdf", path), width = 7.5, height = 10)




p = ggplot(d) +
	facet_grid(mode~type) +
	geom_raster(aes(x = freq - (1/200), y = time, fill = prob), interpolate = T) +
	geom_hline(yintercept = c(10^(0:5)), color = "white", size = 1/2, alpha = 1/3) +
	geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white", size = 1/2, alpha = 1/3) +
	scale_fill_gradientn(colours = rev(c("white", "yellow", "orange", "red", "purple", "blue", "black")), na.value = "grey80", name = "Probability") +
	scale_y_log10(breaks = c(10^(0:5)), labels = c(format(10^(0:5), trim = T, scientific = F, big.mark = ','))) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(expand = F, xlim = c(0,1), ylim = c(0.98, 102000)) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1),
				panel.background = element_rect(fill = "purple"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("TMRCA")
p










for (t in names(gg)) {
	for (f in names(gg[[t]])) {
		if (sum(gg[[t]][[f]]) == 0) {
			gg[[t]][[f]] = NULL
			next
		}
		gg[[t]][[f]] = gg[[t]][[f]] / sum(gg[[t]][[f]])
		gg[[t]][[f]] = as.data.table(as.list(gg[[t]][[f]]))
		gg[[t]][[f]] = cbind(time = as.numeric(t), freq = as.numeric(f), gg[[t]][[f]])
	}
	if (length(gg[[t]]) == 0) {
		gg[[t]] = NULL
		next
	}
	gg[[t]] = rbindlist(gg[[t]])
}
gg = rbindlist(gg)
gg = rbind(data.table(time = gg$time, freq = gg$freq, type = "00", prob = gg$`00`),
					 data.table(time = gg$time, freq = gg$freq, type = "01", prob = gg$`01`),
					 data.table(time = gg$time, freq = gg$freq, type = "02", prob = gg$`02`),
					 data.table(time = gg$time, freq = gg$freq, type = "11", prob = gg$`11`),
					 data.table(time = gg$time, freq = gg$freq, type = "12", prob = gg$`12`),
					 data.table(time = gg$time, freq = gg$freq, type = "22", prob = gg$`22`))
gg = gg[order(gg$time, decreasing = T),]



q = ((0:200)/200)
p = 1 - q

expect.non = rbind(data.table(freq = q, type = "00", prob = p^4),
									 data.table(freq = q, type = "01", prob = 4*(p^3)*q),
									 data.table(freq = q, type = "02", prob = 2*(p^2)*(q^2)),
									 data.table(freq = q, type = "11", prob = 4*(p^2)*(q^2)),
									 data.table(freq = q, type = "12", prob = 4*p*(q^3)),
									 data.table(freq = q, type = "22", prob = q^4))
expect.non$time = 100000
expect.non$mode = "NON"

expect.ibd = rbind(data.table(freq = q, type = "00", prob = p^3),
									 data.table(freq = q, type = "01", prob = 2*(p^2)*q),
									 data.table(freq = q, type = "02", prob = 0),
									 data.table(freq = q, type = "11", prob = ((p^2)*q) + (p*(q^2))),
									 data.table(freq = q, type = "12", prob = 2*p*(q^2)),
									 data.table(freq = q, type = "22", prob = q^3))
expect.ibd$time = 1
expect.ibd$mode = "IBD"

e = rbind(expect.non, expect.ibd)




for (t in names(hh)) {
	for (f in names(hh[[t]])) {
		if (sum(hh[[t]][[f]]) == 0) {
			hh[[t]][[f]] = NULL
			next
		}
		hh[[t]][[f]] = hh[[t]][[f]] / sum(hh[[t]][[f]])
		hh[[t]][[f]] = as.data.table(as.list(hh[[t]][[f]]))
		hh[[t]][[f]] = cbind(time = as.numeric(t), freq = as.numeric(f), hh[[t]][[f]])
	}
	if (length(hh[[t]]) == 0) {
		hh[[t]] = NULL
		next
	}
	hh[[t]] = rbindlist(hh[[t]])
}
hh = rbindlist(hh)
hh = rbind(data.table(time = hh$time, freq = hh$freq, type = "00", prob = hh$`00`),
					 data.table(time = hh$time, freq = hh$freq, type = "01", prob = hh$`01`),
					 data.table(time = hh$time, freq = hh$freq, type = "11", prob = hh$`11`))
hh = hh[order(hh$time, decreasing = T),]



for (t in names(pp)) {
	for (f in names(pp[[t]])) {
		if (sum(pp[[t]][[f]]) == 0) {
			pp[[t]][[f]] = NULL
			next
		}
		pp[[t]][[f]] = pp[[t]][[f]] / sum(pp[[t]][[f]])
		pp[[t]][[f]] = as.data.table(as.list(pp[[t]][[f]]))
		pp[[t]][[f]] = cbind(time = as.numeric(t), freq = as.numeric(f), pp[[t]][[f]])
	}
	if (length(pp[[t]]) == 0) {
		pp[[t]] = NULL
		next
	}
	pp[[t]] = rbindlist(pp[[t]])
}
pp = rbindlist(pp)
pp = rbind(data.table(time = pp$time, freq = pp$freq, type = "00", prob = pp$`00`),
					 data.table(time = pp$time, freq = pp$freq, type = "01", prob = pp$`01`),
					 data.table(time = pp$time, freq = pp$freq, type = "11", prob = pp$`11`))
pp = pp[order(pp$time, decreasing = T),]






library(ggplot2)
library(ggthemes)


p = ggplot(gg) +
	facet_wrap(~type, nrow = 1) +
	geom_point(aes(x = freq - (1/200), y = prob, color = time), size = 2/3, alpha = 2/3, shape = 16) +
	geom_line(data = e, aes(x = freq, y = prob, group = mode), color = "white", size = 4/3, alpha = 1/2) +
	geom_line(data = e, aes(x = freq, y = prob, group = mode), color = "white", size = 2/3) +
	geom_line(data = e, aes(x = freq, y = prob, linetype = mode), size = 2/3) +
	scale_color_gradientn(colours = c("yellow", "orange", "red", "purple", "blue", "black"), trans = "log10", name = "TMRCA",
												breaks = 10^(0:6), labels = format(10^(0:6), trim = T, scientific = F, big.mark = ',')) +
	scale_linetype_manual(values = c("11", "22")) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				panel.background = element_rect(fill = "grey95"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("Probability")
p

ggsave(p, filename = sprintf("_plot.dist.%s.pdf", path), width = 12, height = 4.75)


del = which(gg$time > 100000)
if (length(del) > 0) gg = gg[-del,]

p = ggplot(gg) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(data = e, aes(x = freq, y = time, fill = prob), interpolate = T) +
	geom_raster(aes(x = freq - (1/200), y = time, fill = prob), interpolate = T) +
	geom_hline(yintercept = c(10^(0:5)), color = "white", size = 1/2, alpha = 1/3) +
	geom_hline(yintercept = c(1, 100000), color = "white", size = 2/3) +
	geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white", size = 1/2, alpha = 1/3) +
	scale_fill_gradientn(colours = rev(c("white", "yellow", "orange", "red", "purple", "blue", "black")), na.value = "grey80", name = "Probability") +
	scale_y_log10(breaks = c(10^(0:5)), labels = c(format(10^(0:5), trim = T, scientific = F, big.mark = ','))) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(expand = F, xlim = c(0,1), ylim = c(0.5, 200000)) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1),
				panel.background = element_rect(fill = "grey80"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("TMRCA")
p

ggsave(p, filename = sprintf("_plot.heat.%s.pdf", path), width = 12, height = 4.75)




hh = pp


library(ggplot2)
library(ggthemes)


p = ggplot(hh) +
	facet_wrap(~type, nrow = 1) +
	geom_point(aes(x = freq - (1/200), y = prob, color = time), size = 2/3, alpha = 2/3, shape = 16) +
	scale_color_gradientn(colours = c("yellow", "orange", "red", "purple", "blue", "black"), trans = "log10", name = "TMRCA",
												breaks = 10^(0:6), labels = format(10^(0:6), trim = T, scientific = F, big.mark = ',')) +
	scale_linetype_manual(values = c("11", "22")) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				panel.background = element_rect(fill = "grey95"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("Probability")
p

ggsave(p, filename = sprintf("_plot.H.dist.%s.pdf", path), width = 7.5, height = 5.5)
ggsave(p, filename = sprintf("_plot.P.dist.%s.pdf", path), width = 7.5, height = 5.5)



p = ggplot(hh) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(aes(x = freq - (1/200), y = time, fill = prob), interpolate = T) +
	geom_hline(yintercept = c(10^(0:5)), color = "white", size = 1/2, alpha = 1/3) +
	geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white", size = 1/2, alpha = 1/3) +
	scale_fill_gradientn(colours = rev(c("white", "yellow", "orange", "red", "purple", "blue", "black")), na.value = "grey80", name = "Probability") +
	scale_y_log10(breaks = c(10^(0:5)), labels = c(format(10^(0:5), trim = T, scientific = F, big.mark = ','))) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(expand = F, xlim = c(0,1), ylim = c(0.98, 102000)) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1),
				panel.background = element_rect(fill = "purple"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("TMRCA")
p

ggsave(p, filename = sprintf("_plot.H.heat.%s.pdf", path), width = 7.5, height = 5.5)
ggsave(p, filename = sprintf("_plot.P.heat.%s.pdf", path), width = 7.5, height = 5.5)






p = ggplot(hh) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(aes(x = freq - (1/200), y = time, fill = log10(prob+1e-08)), interpolate = T) +
	geom_hline(yintercept = c(10^(0:5)), color = "white", size = 1/2, alpha = 1/3) +
	geom_hline(yintercept = c(1, 100000), color = "white", size = 2/3) +
	geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white", size = 1/2, alpha = 1/3) +
	scale_fill_gradientn(colours = rev(c("white", "yellow", "orange", "red", "purple", "blue", "black")), na.value = "grey80", name = "log10(P + 1e-8)") +
	scale_y_log10(breaks = c(10^(0:5)), labels = c(format(10^(0:5), trim = T, scientific = F, big.mark = ','))) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(expand = F, xlim = c(0,1), ylim = c(0.5, 200000)) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1),
				panel.background = element_rect(fill = "grey80"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("TMRCA")
p

ggsave(p, filename = sprintf("_plot.logheat.%s.pdf", path), width = 12, height = 4.75)




rr = gg

tru = gg
err = rr

a = sprintf("%0.4f %0.4f %s", tru$time, tru$freq, tru$type)
b = sprintf("%0.4f %0.4f %s", err$time, err$freq, err$type)
i = intersect(a, b)
a = which(a %in% i)
b = which(b %in% i)

tru = tru[a,]
err = err[b,]

tru = tru[order(tru$time, tru$freq, tru$type),]
err = err[order(err$time, err$freq, err$type),]

gg = tru
gg$prob = abs(gg$prob - err$prob)



ggplot(gg) +
	facet_wrap(~type) +
	geom_raster(aes(x = freq - (1/200), y = time, fill = prob), interpolate = T) +
	geom_hline(yintercept = c(10^(0:5)), color = "grey", size = 1/2, alpha = 1/3) +
	geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "grey", size = 1/2, alpha = 1/3) +
	scale_fill_gradientn(colours = (c("white", "red", "black")), na.value = "grey80", name = "Delta", limits = c(0, 0.5)) +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), trim = T, scientific = F, big.mark = ',')) +
	coord_cartesian(expand = F) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1),
				panel.background = element_rect(fill = "grey80"),
				#legend.title = element_blank(),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("TMRCA")








ggplot(gg) +
	facet_wrap(~type) +
	geom_raster(aes(x = time, y = freq, fill = prob+1e-08)) +
	scale_fill_gradientn(colours = (c("yellow", "orange", "red", "purple", "blue", "black")), name = "Probability", trans = "log10") +
	scale_x_log10() 
scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				panel.background = element_rect(fill = "grey95"),
				legend.key.width = unit(2.5, "cm"),
				legend.position = "bottom") +
	xlab("Allele frequency") + ylab("Probability")







cut.time = function(t, brk = c(0, 2*10^(0:5)), lab = sprintf("% 8d", c(10^(0:6)))) {
	as.character(cut(t, breaks = brk, labels = lab, include.lowest = T))
}

cut.geno = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}

make.matrix = function(len = nrow(G)) {
	m = matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
	M = list()
	for (tag in sprintf("% 8d", c(10^(0:6)))) {
		M[[tag]] = m
	}
	M
}


count.pairs = function(ibd, pair, min.size = 10, mask = c("00"=1, "01"=2, "02"=3, "11"=4, "12"=5, "22"=6)) {
	count = make.matrix()
	
	ibd$n = ibd$rhs.index - ibd$lhs.index
	
	del = which(ibd$n < min.size)
	if (length(del) > 0) {
		ibd = ibd[-del, ]
	}
	
	time = cut.time(ibd$time)
	
	for (i in 1:nrow(ibd)) {
		t = time[i]
		
		lhs = ibd$lhs.index[i] + 2
		rhs = ibd$rhs.index[i] - 2
		
		rng = lhs : rhs
		
		coord = matrix(c(rng, mask[pair[rng]]), ncol = 2, byrow = F)
		count[[ t ]][ coord ] = count[[ t ]][ coord ] + 1
	}
	
	count
}


###


cat("Observed genotypes per site ...\n")

count = make.matrix()

i = 0
for (pair in names(ibd)) {
	i = i + 1
	cat(sprintf("%d: %s\n", i, pair))
	
	cmb = ibd[[pair]]
	
	g0 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", pair))
	g1 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", pair))
	g0g1 = cut.geno(G[, g0], G[, g1])
	
	for (tag in names(cmb)) {
		tmp = count.pairs(cmb[[tag]], g0g1)
		for (time in names(tmp)) {
			count[[time]] = count[[time]] + tmp[[time]]
		}
	}
}
cat("\n")


save(count, file = sprintf("result.%s", ibd.file))





###
stop("DONE")
###


library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20"

load(sprintf("../data.%s.RData", prefix))



cut.freq = function(f, brk = (0:101)/101, lab = sprintf("%.3f", (1:101)/101)) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

make.matrix = function(len = length(POS)) {
	m = matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
	M = list()
	for (tag in sprintf("% 8d", c(10^(0:6)))) {
		M[[tag]] = m
	}
	M
}


res.files = dir(pattern = "^result\\..+\\.RData$")

tmp = make.matrix()
for (res.file in res.files) {
	cat(res.file, "\n")
	
	load(res.file)
	
	for (tag in names(count)) {
		tmp[[tag]] = tmp[[tag]] + count[[tag]]
	}
}
count = tmp
rm(tmp)


count$`      10` = count$`      10` + count$`       1`
count$`       1` = NULL

count$`  100000` = count$`  100000` + count$` 1000000`
count$` 1000000` = NULL


for (i in 2:length(count)) {
	count[[i]] = count[[i]] + count[[i-1]]
}


cat("Observed genotypes per frequency bin ...\n")

freq = list()

for (time in names(count)) {
	tmp = as.data.frame(count[[time]])
	tmp = split(tmp, cut.freq(AAF))
	tmp = lapply(tmp, function(x) {
		x = colSums(x)
		x = x / sum(x)
		as.data.table(as.list(x))
	})
	tag = names(tmp)
	tmp = rbindlist(tmp)
	rownames(tmp) = tag
	
	freq[[time]] = tmp
}
cat("\n")



### plotting ...


col = c("00" = "#6082E5",
				"01" = "#46B29D",
				"02" = "#B571E3",
				"11" = "#DBAC12",
				"12" = "#969130",
				"22" = "#DE473A")


d = NULL

for (t in names(freq)) {
	f = as.numeric(rownames(freq[[t]]))
	
	for (p in colnames(freq[[t]])) {
		x = freq[[t]][[p]]
		
		d = rbind(d, data.table(time = t,
														freq = f,
														pair = p,
														prop = x))
	}
}

del = which(is.na(d$prop))
if (length(del) > 0) {
	d = d[-del, ]
}


q = ((0:200)/200)
p = 1 - q

expect.non = rbind(data.table(freq = q, pair = "00", prop = p^4),
									 data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
									 data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
									 data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
									 data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
									 data.table(freq = q, pair = "22", prop = q^4))
expect.non$type = "Expected under NON"

expect.ibd = rbind(data.table(freq = q, pair = "00", prop = p^3),
									 data.table(freq = q, pair = "01", prop = 2*(p^2)*q),
									 data.table(freq = q, pair = "02", prop = 0),
									 data.table(freq = q, pair = "11", prop = ((p^2)*q) + (p*(q^2))),
									 data.table(freq = q, pair = "12", prop = 2*p*(q^2)),
									 data.table(freq = q, pair = "22", prop = q^3))
expect.ibd$type = "Expected under IBD"

e = rbind(expect.non, expect.ibd)
r = cbind(expect.non, ibd = expect.ibd$prop)




d$time[which(d$time == "      10")] = "    t ≤ 10    "
d$time[which(d$time == "     100")] = "   t ≤ 100   "
d$time[which(d$time == "    1000")] = "  t ≤ 1000  "
d$time[which(d$time == "   10000")] = " t ≤ 10000 "
d$time[which(d$time == "  100000")] = "t > 10000"



gg = ggplot(d) +
	facet_wrap(~time, nrow = 1) +
	geom_line(data = expect.ibd, aes(x = freq, y = prop, group = pair), linetype = "33", alpha = 2/3) +
	geom_line(data = expect.non, aes(x = freq, y = prop, group = pair), linetype = "11", alpha = 2/3) +
	geom_ribbon(data = r, aes(x=freq, ymin = prop, ymax = ibd, fill = pair), alpha = 1/3) +
	geom_point(aes(x = freq, y = prop, colour = pair), alpha = 0.8, size = 1) +
	coord_cartesian(xlim = c(-0.025, 1.025), ylim = c(-0.01, 1.01), expand = F) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_colour_manual(values = col) +
	scale_fill_manual(values = col) +
	theme_few() +
	theme(aspect.ratio = 2.8125,
				panel.border = element_rect(fill = NA, colour = "black", size = 1/2),
				legend.title = element_blank(),
				legend.background = element_rect(fill = "white", colour = "black", size = 1/3),
				legend.position = c(1 - (1/10.5),0.825),
				panel.grid = element_line(colour = "grey90", size = 1/2),
				panel.grid.major = element_line(colour = "grey90", size = 1/2),
				panel.grid.minor = element_blank()) +
	xlab("Allele frequency") + ylab("Genotype pair proportion") +
	guides(color = guide_legend(override.aes = list(alpha = 1)))
gg

ggsave(filename = "./_plot.genpair_tmrca.pdf", plot = gg, width = 11, height = 6.5)


ggplot(d) +
	facet_wrap(~pair) +
	geom_point(aes(x = freq, y = prop, colour = time), size = 0.25)













