#
# time & freq dependend genotype pair emissions
#

library(data.table)

args = commandArgs(T)
prefix = args[1]

cat("Loading...")

load(sprintf("./data.%s.RData", prefix))
#load(sprintf("./data.%s.G.RData", prefix)) # genotypes
load(sprintf("./data.%s.H.RData", prefix)) # haplotypes
#load(sprintf("./data.%s.P.RData", prefix)) # phased haplotypes

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
combine.rng = function(mat, len = nrow(H)) {
	#combine.rng = function(mat, len = nrow(P)) {
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
	
	save.file = sprintf("result_P.%s", basename(ibd.file))
	flag.file = sprintf("_flag_P.%s.txt", basename(ibd.file))
	
	if (file.exists(save.file)) next
	if (file.exists(flag.file)) next
	
	cat(".", file = flag.file)
	
	#GG = make.gen.array()
	#HH = make.hap.array()
	PP = make.hap.array()
	
	load(ibd.file)
	
	i = 0
	for (pair in names(ibd)) {
		i = i + 1
		cat(sprintf("%d: %s\n", i, pair))
		
		cmb = ibd[[pair]]
		
		# 		g0 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", pair))
		# 		g1 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", pair))
		#
		# 		a = G[, g0]
		# 		b = G[, g1]
		# 		gg = sprintf("%d%d", a, b)
		# 		x = which(gg == "10"); if (length(x) > 0) gg[x] = "01"
		# 		x = which(gg == "20"); if (length(x) > 0) gg[x] = "02"
		# 		x = which(gg == "21"); if (length(x) > 0) gg[x] = "12"
		
		
		cmb = lapply(cmb, prune.ibd)
		
		for (haps in names(cmb)) {
			mat = cmb[[haps]]
			
			h0 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", haps))
			h1 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", haps))
			
			# a = H[, h0]
			# b = H[, h1]
			# hh = sprintf("%d%d", a, b)
			# x = which(hh == "10"); if (length(x) > 0) hh[x] = "01"
			
			a = P[, h0]
			b = P[, h1]
			pp = sprintf("%d%d", a, b)
			x = which(pp == "10"); if (length(x) > 0) pp[x] = "01"
			
			crd = cut.time(mat$time)
			crd = split(mat, crd)
			crd = lapply(crd, combine.rng)
			crd = lapply(crd, function(x, f) {
				#split((1:nrow(G))[x], f[x])
				#split((1:nrow(H))[x], f[x])
				split((1:nrow(P))[x], f[x])
			}, FRQ)
			
			for (t in names(crd)) {
				for (f in names(crd[[t]])) {
					rng = crd[[t]][[f]]
					
					if (length(rng) == 0) next
					
					#g = table(gg[rng])
					#h = table(hh[rng])
					p = table(pp[rng])
					
					# for (tag in names(g)) {
					# 	GG[[t]][[f]][tag] = GG[[t]][[f]][tag] + g[tag]
					# }
					# for (tag in names(h)) {
					# 	HH[[t]][[f]][tag] = HH[[t]][[f]][tag] + h[tag]
					# }
					for (tag in names(p)) {
						PP[[t]][[f]][tag] = PP[[t]][[f]][tag] + p[tag]
					}
				}
			}
		}
	}
	
	#save(GG, HH, PP, file = save.file)
	save(PP, file = save.file)
	unlink(flag.file)
}





stop()

library(data.table)
library(ggplot2)
library(ggthemes)


expected.age = function(x) {
	x0 = which(x == 0)
	x1 = which(x == 1)
	ex = ((-2 * x) / (1 - x)) * log(x)
	if (length(x0) > 0) ex[x0] = 0
	if (length(x1) > 0) ex[x1] = 2
	ex
}


path = "tru"
path = "err"

#gg = make.gen.array()
hh = make.hap.array()
#pp = make.hap.array()

res.files = dir(pattern = "^result_H\\..+", path = path, full.names = T)

for (res.file in res.files) {
	print(res.file)
	load(res.file)
	
	# for (t in names(GG)) {
	# 	for (f in names(GG[[t]])) {
	# 		gg[[t]][[f]] = gg[[t]][[f]] + GG[[t]][[f]]
	# 	}
	# }
	for (t in names(HH)) {
		for (f in names(HH[[t]])) {
			hh[[t]][[f]] = hh[[t]][[f]] + HH[[t]][[f]]
		}
	}
	# for (t in names(PP)) {
	# 	for (f in names(PP[[t]])) {
	# 		pp[[t]][[f]] = pp[[t]][[f]] + PP[[t]][[f]]
	# 	}
	# }
}

#hh = pp

# hh = lapply(hh, function(x) {
# 	lapply(x, function(y) {
# 		y / sum(y)
# 	})
# })

se <- function(x) sqrt(var(x, na.rm = T)/length(which(!is.na(x))))


dd = NULL


do.ibd = F
d = NULL

do.ibd = T
d = NULL


xx = lapply(hh, function(x) as.matrix(t(as.data.frame(x))))
ts = as.numeric(names(hh))
#ts = ts[1:80]
#ks = c(expected.age( exp(seq(log(0.001), log(1), length.out = 11)) ) * 7300, 4*7300 ) # , 100000) #, 1000000)
ne = 7300
ks = c(0, 1/4, 1/2, 1, 2, 4, 8) * ne
ks = c(0, 100, 1000, 1000)
for (k in 2:(length(ks)-1)) {
	t = NULL
	
	if (do.ibd) {
		#t = which(ts > ks[k-1] & ts <= ks[k])
		#t = which(ts > ks[1] & ts <= ks[k])
		t = which(ts <= ks[k])
		#t = which(ts >= ks[k-1] & ts <= ks[k+1])
		#t = which.min(abs(ts - ks[k]))
		#t = c(t-1, t, t+1)
	} else {
		t = which(ts > ks[k])
		#t = which(ts >= ks[k] & ts < ks[k+1])
		#t = 1:length(ts)
	}
	
	# g0 = lapply(xx[t], function(x) data.table(x[, 1]  ))
	# g1 = lapply(xx[t], function(x) data.table(x[, 2]  ))
	# g2 = lapply(xx[t], function(x) data.table(x[, 3]  ))
	# 
	# g0 = as.matrix(Reduce(cbind, g0))
	# g1 = as.matrix(Reduce(cbind, g1))
	# g2 = as.matrix(Reduce(cbind, g2))
	# 
	# f0 = g0 / (g0+g1+g2)
	# f1 = g1 / (g0+g1+g2)
	# f2 = g2 / (g0+g1+g2)
	# 
	# mn = data.table("00" = apply(f0, 1, mean, na.rm = T), "01" = apply(f1, 1, mean, na.rm = T), "11" = apply(f2, 1, mean, na.rm = T))
	# er = data.table("00" = apply(f0, 1, se), "01" = apply(f1, 1, se), "11" = apply(f2, 1, se))
	
	x = Reduce("+", xx[t])
	mn = t(apply(x, 1, function(q) q / sum(q)))
	x = as.data.table(mn)
	
	tmp = NULL
	for (i in 1:ncol(mn)) {
		tm = NULL
		if (do.ibd) {
			tm = paste("≤", ks[k])
			#tmp = rbind(tmp, data.table(time = sprintf("%s - %s N[e]", format(round(ks[k-1]/ne, 2), trim = T), format(round(ks[k]/ne, 2), trim = T)), freq = (seq(0, 1, length.out = 101))[-1] - 0.005, type = names(x)[i], prob = x[[i]]))
		} else {
			tm = paste(">", ks[k])
			#tmp = rbind(tmp, data.table(time = sprintf("%s - %s N[e]", format(round(ks[k]/ne, 2), trim = T), format(round(ks[k+1]/ne, 2), trim = T)), freq = (seq(0, 1, length.out = 101))[-1] - 0.005, type = names(x)[i], prob = x[[i]]))
		}
		
		# tmp = rbind(tmp, data.table(time = ks[k], freq = (seq(0, 1, length.out = 101))[-1] - 0.005, type = names(mn)[i], prob = data.frame(mn)[,i], err = data.frame(er)[,i]))
		tmp = rbind(tmp, data.table(time = ks[k], freq = (seq(0, 1, length.out = 101))[-1] - 0.005, type = names(x)[i], prob = x[[i]]))
		
	}
	d = rbind(d, tmp)
}
del = which(is.na(d$prob))
if (length(del) > 0) d = d[-del,]

# include expected
d = rbind(data.table(time = unique(d$time), freq = 0, type = "00", prob = 1),
					data.table(time = unique(d$time), freq = 0, type = "01", prob = 0),
					data.table(time = unique(d$time), freq = 0, type = "11", prob = 0),
					data.table(time = unique(d$time), freq = 1, type = "00", prob = 0),
					data.table(time = unique(d$time), freq = 1, type = "01", prob = 0),
					data.table(time = unique(d$time), freq = 1, type = "11", prob = 1),
					d)
d = as.data.frame(d)
d = unique(d)


if (do.ibd)  d$mode = "< TMRCA"
if (!do.ibd) d$mode = "> TMRCA"

dd = rbind(dd, d)




dd = dd[order(dd$mode, dd$time, dd$freq, dd$type),]


tmp = split(dd, list(dd$time, dd$type, dd$mode))
tmp = lapply(tmp, function(d) {
	#w = c(1e16, 1e8, rep(1, nrow(d)-4), 1e8, 1e16)
	l = loess(prob ~ freq, data = d) # , weights = w)
	l = predict(l, data.frame(freq = d$freq), se = TRUE)
	d$fit = l$fit
	x = which(d$fit < 0); if (length(x) > 0) d$fit[x] = 0
	x = which(d$fit > 1); if (length(x) > 0) d$fit[x] = 1
	d
})
ll = rbindlist(tmp)

tmp = split(ll, list(ll$time, ll$mode))
tmp = lapply(tmp, function(d) {
	x = split(d, d$type)
	sum = x[[1]]$fit + x[[2]]$fit + x[[3]]$fit
	x[[1]]$fit = x[[1]]$fit /sum
	x[[2]]$fit = x[[2]]$fit /sum
	x[[3]]$fit = x[[3]]$fit /sum
	rbindlist(x)
})
ll = rbindlist(tmp)


ll$outl = F
x = which(abs(ll$fit - ll$prob) > 0.025)
if (length(x) > 0) {
	y = which(ll$freq[x] == 0 | ll$freq[x] == 1)
	if (length(y) > 0) x = x[-y]
	ll$outl[x] = T
}


# tmp = split(ll, list(ll$time, ll$type, ll$mode))
# tmp = lapply(tmp, function(z) {
# 	f = z$freq
# 	p = z$prob
# 	n = nrow(z)
# 	q = p
# 	for (j in 1:100) {
# 		for (i in 3:(n-2)) {
# 			q[i] = mean(c(p[i-1], p[i], p[i+1]))
# 		}
# 		p = q
# 	}
# 	z$prob = p
# 	z
# })
# lll = rbindlist(tmp)

# tmp = split(lll, list(lll$time, lll$mode))
# tmp = lapply(tmp, function(d) {
# 	x = split(d, d$type)
# 	sum = x[[1]]$fit + x[[2]]$fit + x[[3]]$fit
# 	x[[1]]$fit = x[[1]]$fit /sum
# 	x[[2]]$fit = x[[2]]$fit /sum
# 	x[[3]]$fit = x[[3]]$fit /sum
# 	rbindlist(x)
# })
# lll = rbindlist(tmp)



p = ggplot(ll) +
	facet_grid(mode~type) +
	#geom_point(aes(x = freq, y = prob, color = factor(time), shape = outl), size = 1.5, alpha = 2/3) +
	geom_line(aes(x = freq, y = fit, color = factor(time)), alpha = 2/3, show.legend = F) +
	geom_point(aes(x = freq, y = prob, color = factor(time)), size = 1/3, alpha = 1, show.legend = F) +
	#geom_linerange(aes(x = freq, ymin = prob-err, ymax = prob+err, color = factor(time)), alpha = 1/2, show.legend = F) +
	#geom_smooth(data = ee, aes(x = freq, y = prob, color = time), method = "loess", size = 1/2) + # , method.args = list(family = "symmetric", surface="interpolate", iterations=1, cell = 1)) +
	#scale_color_manual(values = c("yellow", "orange", "red", "purple", "blue", "black"), name = "TMRCA") +
	#scale_color_manual(values = c("orange", "red", "blue"), name = "TMRCA") +
	#scale_linetype_manual(values = c("11", "22")) +
	scale_shape_manual(values = c(16, 4), guide=FALSE) +
	#scale_color_manual(values = c("darkorange", "limegreen", "dodgerblue", "purple1")) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
	coord_cartesian(xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				strip.text = element_text(face = "bold", size = 11),
				plot.title = element_text(face = "bold"),
				axis.text = element_text(size = 8),
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				panel.background = element_rect(fill = "grey95")) +
				#legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2)) +
				#legend.direction = "horizontal",
				#legend.position = "bottom") + # c(0.5,0.845)) +
	xlab("Allele frequency") + ylab("Probability") +
	guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
	ggtitle("Before error")
p


save(ll, file = sprintf("_data.H.emiss.%s.RData", path))
ggsave(p, filename = sprintf("_plot.H.emiss.%s.pdf", path), width = 7, height = 8.5)


ggsave(p, filename = sprintf("__plot.H.exp.%s.pdf", path), width = 4.5*(1), height = 5.5*(1))

ggsave(p, filename = sprintf("__plot.H.exp.%s.exl.pdf", path), width = 4.5*(1), height = 5.5*(1))


ggsave(p, filename = sprintf("__plot.P.exp.%s.pdf", path), width = 7.5, height = 5.5)




path = "err"
load(sprintf("_data.H.emiss.%s.RData", path))


tmp = split(ll, list(ll$time, ll$type, ll$mode))
tmp = lapply(tmp, function(d) {
	f = (0:500)/500
	ye = approx(d$freq, d$prob, xout = f)$y
	yf = approx(d$freq, d$fit,  xout = f)$y
	data.table(time = d$time[1], freq = f, type = d$type[1], prob = ye, fit = yf, mode = d$mode[1])
})
z = rbindlist(tmp)
z = split(z, list(z$time, z$mode))
z = lapply(z, function(d) {
	x = split(d, d$type)
	sum = x[[1]]$fit + x[[2]]$fit + x[[3]]$fit
	x[[1]]$fit = x[[1]]$fit /sum
	x[[2]]$fit = x[[2]]$fit /sum
	x[[3]]$fit = x[[3]]$fit /sum
	sum = x[[1]]$prob + x[[2]]$prob + x[[3]]$prob
	x[[1]]$prob = x[[1]]$prob /sum
	x[[2]]$prob = x[[2]]$prob /sum
	x[[3]]$prob = x[[3]]$prob /sum
	rbindlist(x)
})
z = rbindlist(tmp)


for (t in sort(unique(z$time))) {
	x = z[which(z$time == t), ]
	
	non = which(x$mode == "> TMRCA")
	non = x[non,]
	ibd = which(x$mode == "< TMRCA")
	ibd = x[ibd,]
	
	if (!identical(non$freq, ibd$freq)) stop("???")
	
	a = sapply(split(non, non$type), function(x) x$prob)
	colnames(a) = sprintf("NON_%s", colnames(a))
	a = apply(a, 2, function(x) round(x/rowSums(a), 8))
	a = apply(a, 2, function(x) format(x, scientific = F))
	
	b = sapply(split(ibd, ibd$type), function(x) x$prob)
	colnames(b) = sprintf("IBD_%s", colnames(b))
	b = apply(b, 2, function(x) round(x/rowSums(b), 8))
	b = apply(b, 2, function(x) format(x, scientific = F))
	
	out = cbind(Frequency = format(unique(x$freq), scientific = F), a, b)
	
	write.table(out, file = sprintf("_emission.hhmm.%s.%s.txt", path, t), append = F, sep = " ", quote = F, row.names = F, col.names = T)
}




stop()


### initials


load("~/Research/DetectIBD/data.OutOfAfricaHapMap20.H.RData")
tru = H

load("~/Research/DetectIBD/GenErr_1000G/data.OutOfAfricaHapMap20.GenErr_1000G.H.RData")
err = H

n = nrow(err)
s = ncol(err)
ss = 1:s

ft = rowSums(tru)
fe = rowSums(err)

ii = 1:nrow(err)
# ii = split(ii, fe)
# ii = lapply(ii, function(x) {
# 	if (length(x) < 1000) return(x)
# 	sample(x, 1000)
# })
# ii = sort(as.vector(unlist(ii)))
# ff = fe[ii]
# nn = length(ii)

con = matrix(NA, nrow = n, ncol = 3)
dis = matrix(NA, nrow = n, ncol = 3)


open = rep(F, n)
open[which(is.na(con[, 1]) | is.na(dis[, 1]))] = T
j = 0
for (i in sample(ii)) {
	if (fe[i] < 4000 | fe[i] > s-2) next
	if (!open[i]) next
	j = j + 1
	cat(sprintf("(%d) %d of %d\n", j, i, n))
	
	x = which(err[i, ] == 1)
	
	cmb = combn(x, 2)
	c1 = as.matrix(data.frame(i, cmb[1,]))
	c2 = as.matrix(data.frame(i, cmb[2,]))
	
	d1 = as.matrix(data.frame(i, sample(x, s, replace = T)))
	d2 = as.matrix(data.frame(i, sample(setdiff(ss, x), s, replace = T)))
	
	con[i,1] = fe[i]
	con[i,2] = length(which(tru[ c1 ] + tru[ c2 ] == 2))
	con[i,3] = ncol(cmb)
	
	dis[i,1] = fe[i]
	dis[i,2] = length(which(tru[ d1 ] + tru[ d2 ] == 1))
	dis[i,3] = s
}


save(con, dis, file = "_initial.sampling.RData")
load("_initial.sampling.RData")


#load("~/Dropbox/tmp/hhmm/__tmp.RData")
se <- function(x) sqrt(var(x)/length(x))

C = data.table(f = con[, 1], x = con[, 2], n = con[, 3])
C = C[-(which(is.na(C$x))),]
C$p = C$x / C$n

cc = split(C, C$f)
cc = lapply(cc, function(x) {
	data.table(f = x$f[1], m = mean(x$p), s = se(x$p), n = nrow(x))
})
cc = rbindlist(cc)

#cc = cc[-(which(cc$n < 2)),]
plot(cc$f, cc$m)


D = data.table(f = dis[, 1], x = dis[, 2], n = dis[, 3])
D = D[-(which(is.na(D$x))),]
D$p = D$x / D$n

dd = split(D, D$f)
dd = lapply(dd, function(x) {
	data.table(f = x$f[1], m = mean(x$p), s = se(x$p), n = nrow(x))
})
dd = rbindlist(dd)

plot(dd$f, dd$m)


ggplot(dd) + 
	#geom_pointrange(aes(x=f/5000, y=m, ymin=m-s, ymax=m+s), size = 1/10, alpha = 1/3) +
	geom_pointrange(aes(x=f, y=m, ymin=m-s, ymax=m+s), size = 1/10, alpha = 1/3) +
	#coord_cartesian(xlim=c(-0.01, 1.01), ylim=c(0.799, 1.001), expand = F) +
	coord_cartesian(xlim=c(1.1, 50.9), ylim=c(0.799, 1.001), expand = F) +
	scale_x_continuous(breaks=c(2, 5, 10, 20, 30, 40, 50)) +
	theme_bw() +
	xlab("Allele count, k") + 
	ylab("Empirical initial probability, mean (±SE)") +
  #ggtitle("Concordant pairs")
	ggtitle("Discordant pairs")
	
ggsave("_plot.inits.con.pdf", width = 8, height = 4.5)
ggsave("_plot.inits.dis.pdf", width = 8, height = 4.5)
ggsave("_plot.inits.con.zoom.pdf", width = 3, height = 4.5)
ggsave("_plot.inits.dis.zoom.pdf", width = 3, height = 4.5)


x = which(cc$f < 2 | cc$f > s-2)
if (length(x) > 0) cc = cc[-x,]

x = which(dd$f < 2 | dd$f > s-2)
if (length(x) > 0) dd = dd[-x,]

x = which(is.na(cc$m))
if (length(x) > 0) cc = cc[-x,]

x = which(is.na(dd$m))
if (length(x) > 0) dd = dd[-x,]


k = sort(intersect(cc$f, dd$f))
cc = cc[sort(which(cc$f %in% k)), ]
dd = dd[sort(which(dd$f %in% k)), ]

identical(cc$f, dd$f)

out = data.table(Frequency = format(round(cc$f / s, 4), scientific = F),
								 CON_NON = format(round(1 - cc$m, 8), digits=8, scientific = F),
								 CON_IBD = format(round(    cc$m, 8), digits=8, scientific = F),
								 DIS_NON = format(round(1 - dd$m, 8), digits=8, scientific = F),
								 DIS_IBD = format(round(    dd$m, 8), digits=8, scientific = F))

write.table(out, file = sprintf("_initials.hhmm.%s.txt", path), append = F, sep = " ", quote = F, row.names = F, col.names = T)


out$CON_NON = format(1e-08, digits=8, scientific = F)
out$CON_IBD = format(1-1e-08, digits=8, scientific = F)
out$DIS_NON = format(1e-08, digits=8, scientific = F)
out$DIS_IBD = format(1-1e-08, digits=8, scientific = F)

write.table(out, file = sprintf("_initials.hhmm.%s.txt", "tru"), append = F, sep = " ", quote = F, row.names = F, col.names = T)












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







B = function(a) {
	prod(gamma(a+1)) / gamma(sum(a+1))
}

D = function(x, a) {
	prod(x^a) / B(a)
}




