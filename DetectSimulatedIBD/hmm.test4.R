#
# run HMM tests
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)

###

expected.age = function(k, n) {
	if (length(k) > 1) return(sapply(k, expected.age, n))
	j = 2:n
	s = sum( choose(n - j, k - 1) * ((n - j + 1) / (n * (j - 1))) )
	2 * choose(n - 1, k)^-1 * s
}

expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}


# get.init.prob = function(k, n) {
# 	t = expected.age(k, n)
# 	p = (1/2) * (1 - exp(-6 * t))
# 	#p = 1e-5
# 	c(IBD = p, NON = 1 - p)
# }

get.init.prob = function(k, empirical.init) {
	#p = empirical.init$ibd[ which(empirical.init$fk == k) ]
	p = 1 - 1e-5
	c(IBD = p, NON = 1 - p)
}

get.trans.prob = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
# 	t = expected.age(k, n)
# 	r = d * rec.rate
# 	#p = (1/2) * (1 - exp(-6 * t))
# 	
# 	q11 = exp(-4 * N * r * t)
# 	q10 = 1 - q11
# 	q01 = 0 # (p / (1 - p)) * q10
# 	q00 = 1 # 1 - q01
	
	p = 1e-5
	
	q10 = 5e-8
	q11 = 1 - q10
	q01 = (p / (1 - p)) * q10
	q00 = 1 - q01
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}

# plot trans prob
g = NULL
n = 5000
for (d in c(0, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000)) {
	for (k in c(2, 5, 10, 15, 20, 25)) {
		tmp = as.data.table(melt(get.trans.prob(d, k, n)))
		names(tmp) = c("from", "to", "prob")
		tmp$d = d
		tmp$k = k
		g = rbind(g, tmp)
	}
}

ggplot(g) + 
	facet_grid(to~from) + 
	geom_line(aes(x=d, y=prob, colour = factor(k))) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Inter-variant distance, d, in basepair units") + 
	ylab("Transition probability")


get.emiss.prob = function(k, empirical.ibd, empirical.non) {
	matrix(c(empirical.ibd[k + 1, ], empirical.non[k + 1, ]), nrow = 2, byrow = T, dimnames = list(c("IBD", "NON"), colnames(empirical.ibd)))
}


###


hmm.viterbi = function (obs, init, emiss, trans, cutoff = NULL, return.all = F) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	v = array(NA, c(n.sts, n.obs))
	w = array(NA, n.obs)
	z = array(NA, n.obs)
	
	for (state in states) {
		v[state, 1] = init[state] * emiss[state, obs[1], 1]
	}
	z[1] = 1 #emiss[1, obs[1], 1]
	
	w[1] = max(v[, 1])
	v[, 1] = v[, 1] / w[1]
	#z[1] = z[1] / w[1]
	
	for (i in 2:n.obs) {
		o = obs[i]
		for (state in states) {
			max = 0
			for (prev.state in states) {
				tmp = v[prev.state, i - 1] * trans[prev.state, state, i]
				max = max(max, tmp)
			}
			v[state, i] = emiss[state, o, i] * max
		}
		z[i] = z[i - 1] * trans[1, 1, i] * emiss[1, o, i]
		
		w[i] = max(v[, i])
		v[, i] = v[, i] / w[i]
		z[i] = z[i] / w[i]
		
		if (!is.null(cutoff) && (z[i] / v[1, i]) < cutoff) {
			n.obs = i
			break
		}
	}
	
	path = rep(NA, n.obs)
	
	tmp = which.max(v[, n.obs])
	path[i] = states[tmp]
	
	for (i in (n.obs - 1):1) {
		tmp = which.max(v[ , i] * trans[ , path[i + 1], i])
		path[i] = states[tmp]
	}
	
	if (return.all) {
		return(list(path = path, weights = w, halt = z))
	}
	
	return(path)
}


hmm.fwd = function (obs, init, emiss, trans) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	f = array(NA, c(n.sts, n.obs))
	w = array(NA, n.obs)
	
	for (state in states) {
		f[state, 1] = init[state] * emiss[state, obs[1], 1]
	}
	
	w[1] = max(f[, 1])
	f[, 1] = f[, 1] / w[1]
	
	for (i in 2:n.obs) {
		for (state in states) {
			tmp = 0
			for (prev.state in states) {
				tmp = tmp + (f[prev.state, i - 1] * trans[prev.state, state, i])
			}
			f[state, i] = emiss[state, obs[i], i] * tmp
		}
		
		w[i] = max(f[, i])
		f[, i] = f[, i] / w[i]
	}
	
	return(list(f=f, w=w))
}	


hmm.bwd = function (obs, init, emiss, trans) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	b = array(NA, c(n.sts, n.obs))
	w = array(NA, n.obs)
	
	for (state in states) {
		b[state, n.obs] = 1
	}
	
	w[n.obs] = 1
	
	for (i in (n.obs - 1):1) {
		for (state in states) {
			tmp = 0
			for (next.state in states) {
				tmp = tmp + (b[next.state, i + 1] * trans[state, next.state, i + 1] * emiss[next.state, obs[i + 1], i + 1])
			}
			b[state, i] = tmp
		}
		
		w[i] = max(b[, i])
		b[, i] = b[, i] / w[i]
	}
	
	return(list(b=b, w=w))
}


hmm.posterior = function (obs, init, emiss, trans) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	fwd = hmm.fwd(obs, init, emiss, trans)
	bwd = hmm.bwd(obs, init, emiss, trans)
	
	logl.state = array(NA, c(n.sts, n.obs))
	#logl.model = array(0, n.obs)
	
	wght.fwp = cumsum(log(fwd$w))
	wght.bwp = rev(cumsum(rev(log(bwd$w))))
	
	for (state in states) {
		logl.state[state, ] = log(fwd$f[state, ]) + log(bwd$b[state, ]) + wght.fwp + wght.bwp
		#logl.model = logl.model + (fwd$f[state, ] * bwd$b[state, ])
	}
	
	#logl.model = log(logl.model) + wght.fwp + wght.bwp
	
	logp.fwd = wght.fwp[n.obs] + log(fwd$f[1, n.obs])
	#logp.bwd = wght.bwp[  1  ] + log(bwd$b[1,   1  ] * emiss[1, obs[  1  ],   1  ] * init[1])
	
	for (i in 2:n.sts) {
		logp.fwd = logp.fwd + log(1 + exp((wght.fwp[n.obs] + log(fwd$f[i, n.obs])) - logp.fwd))
		#logp.bwd = logp.bwd + log(1 + exp((wght.bwp[  1  ] + log(bwd$b[i,   1  ] * emiss[i, obs[  1  ],   1  ] * init[i])) - logp.bwd))
	}
	
	return(exp(logl.state - as.numeric(logp.fwd)))
}




###


prefix = "history.generror_1000g"


load(sprintf("hmm.probs.%s.RData", prefix))

load(sprintf("%s.RData", prefix))

load("../truth/ibd_mrca.truth.RData")


sample.size = 5000
effect.size = 10000
recomb.rate = 1e-8

# genotype data

G = as.matrix(load.bigmatrix(sprintf("%s.G", prefix)))

AC = rowSums(G)
AF = AC / (ncol(G) * 2)


# truth data

# del = union(which(truth$wall.lhs), which(truth$wall.rhs))
# if (length(del) > 0) {
# 	truth = truth[-del, ]
# }

del = which(FKI$frq > 0.5)
if (length(del) > 0) {
	del = FKI$index[del]
	del = which(truth$index %in% del)
	if (length(del) > 0) {
		truth = truth[-del, ]
	}
}


# preparing emission probs
# emiss = array(NA, c(2, 6, length(POS)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
# for (i in 1:length(POS)) {
# 	emiss[1 , , i] = hmm.emiss.ibd[AC[i] + 1, ]
# 	emiss[2 , , i] = hmm.emiss.non[AC[i] + 1, ]
# }
emiss.all = array(NA, c(2, 6, length(POS)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
for (i in 1:length(POS)) {
	emiss.all[1 , , i] = hmm.emiss.ibd[AC[i] + 1, ] #+ 1e-8
	emiss.all[2 , , i] = hmm.emiss.non[AC[i] + 1, ] #+ 1e-8
}


rare.vars = which(AC >=2 & AC <= 25)


# random test

fki = FKI[sample(which(FKI$n.sharer == 3), 1), ]

g.sharers = as.numeric(strsplit(fki$g.sharer, "|", T)[[1]])
g.sharers = combn(g.sharers, 2)

h.sharers = as.numeric(strsplit(fki$h.sharer, "|", T)[[1]])
h.sharers = combn(h.sharers, 2)


#init = get.init.prob(AC[fki$index], sample.size)
init = get.init.prob(AC[fki$index], hmm.init)

dist.lhs = c(1, diff(POS))
dist.rhs = c(diff(POS), 1)

trans.lhs = sapply(1:length(POS), function(x) get.trans.prob(dist.lhs[x], AC[fki$index], sample.size), simplify = "array")
trans.rhs = sapply(1:length(POS), function(x) get.trans.prob(dist.rhs[x], AC[fki$index], sample.size), simplify = "array")



d.brks = NULL
d.true = NULL
d.rare = NULL
d.othr = NULL
d.path = NULL
d.post = NULL
d.focl = NULL
d.emis = NULL
d.halt = NULL

for (share in 1:ncol(g.sharers)) {
	
	pair.g = g.sharers[, share]
	pair.h = h.sharers[, share]
	
	g0 = G[, pair.g[1]]
	g1 = G[, pair.g[2]]
	
	foc = fki$index
	min = 1
	max = length(POS)
	
	cat("fk =", AC[foc], "\n")
	
	if (g0[foc] != 1 || g0[foc] != 1) {
		cat("Rare variant not detected!\n")
		next
	}
	
	tag = paste(pair.g, collapse = " + ")
	
	
	# fetch truth
	
	tru = which(truth$index == foc & truth$h0 == pair.h[1] & truth$h1 == pair.h[2])
	
	if (length(tru) == 1) {
		true.position = c(truth$lhs[tru], truth$rhs[tru])
		true.index = match(true.position, POS)
		
		d.true = rbind(d.true, data.table(tag = tag,  
																			idx = true.index,  
																			pos = true.position))
		
		cat("True LHS rel. index:", foc - true.index[1], "\n")
		cat("True RHS rel. index:", true.index[2] - foc, "\n")
		
		rare = truth$index[ which(truth$h0 == truth$h0[tru] & truth$h1 == truth$h1[tru]) ]
		
		d.rare = rbind(d.rare, data.table(tag = tag,  
																			idx = rare,  
																			pos = POS[rare]))
	}
	
	
	# observed genotye pairs
	
	obs.all = sprintf("%d%d", g0, g1)
	i = which(obs.all == "10"); if (length(i) > 0) obs.all[i] = "01"
	i = which(obs.all == "20"); if (length(i) > 0) obs.all[i] = "02"
	i = which(obs.all == "21"); if (length(i) > 0) obs.all[i] = "12"
	
	lhs.rng = foc : min
	rhs.rng = foc : max
	
	lhs.obs = obs.all[ lhs.rng ]
	rhs.obs = obs.all[ rhs.rng ]
	
	
	# emission profile
	e.ibd = rep(NA, length(obs.all))
	e.non = rep(NA, length(obs.all))
	for (i in 1:length(obs.all)) {
		e.ibd[i] = emiss.all["IBD", obs.all[i], i]
		e.non[i] = emiss.all["NON", obs.all[i], i]
	}
	d.emis = rbind(d.emis, data.table(tag = tag, idx = 1:length(POS), pos = POS, ibd = e.ibd, non = e.non))
	
	
	# emiss probs
	lhs.emiss = emiss.all[ , , lhs.rng ]
	rhs.emiss = emiss.all[ , , rhs.rng ]
	
	# trans probs
	lhs.trans = trans.lhs[ , , lhs.rng ]
	rhs.trans = trans.rhs[ , , rhs.rng ]
	
	
	# HMM: Viterbi path
	cat("Running full Viterbi ... \n")
	lhs.path = hmm.viterbi(lhs.obs, init, lhs.emiss, lhs.trans)
	rhs.path = hmm.viterbi(rhs.obs, init, rhs.emiss, rhs.trans)
	
	path = c(rev(lhs.path[-1]), rhs.path)
	i = which(path == 1); if (length(i) > 0) path[i] = "IBD"
	i = which(path == 2); if (length(i) > 0) path[i] = "NON"
	
	d.path = rbind(d.path, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,  state = path))
	
	
	cat("Running short Viterbi ... \n")
	lhs.halt = hmm.viterbi(lhs.obs, init, lhs.emiss, lhs.trans, cutoff = 1e-16)
	rhs.halt = hmm.viterbi(rhs.obs, init, rhs.emiss, rhs.trans, cutoff = 1e-16)
	
	halt = c(rev(lhs.halt[-1]), rhs.halt)
	i = which(halt == 1); if (length(i) > 0) halt[i] = "IBD"
	i = which(halt == 2); if (length(i) > 0) halt[i] = "NON"
	tmp = (foc - length(lhs.halt) + 1) : (foc + length(rhs.halt) - 1)
	d.halt = rbind(d.halt, data.table(tag = tag,  idx = tmp,  pos = POS[tmp],  state = halt))
	
	
	# HMM: Posteriors
	cat("Running posteriors ... \n")
	lhs.post = hmm.posterior(lhs.obs, init, lhs.emiss, lhs.trans)
	rhs.post = hmm.posterior(rhs.obs, init, rhs.emiss, rhs.trans)

	post = array(NA, c(2, length(obs.all)))
	post[1, ] = c(rev(lhs.post[1, -1]), rhs.post[1, ])
	post[2, ] = c(rev(lhs.post[2, -1]), rhs.post[2, ])
	
	d.post = rbind(d.post, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,  pprob = post[1, ]))
	
	
	brks = which(obs.all == "02")
	d.brks = rbind(d.brks, data.table(tag = tag,  idx = brks,  pos = POS[brks]))
	
	d.focl = rbind(d.focl, data.table(tag = tag,  idx = foc,  pos = POS[foc]))
	
	# other rare vars shared
	is.shared = which(g0[rare.vars] == 1 & g1[rare.vars] == 1)
	rare = rare.vars[is.shared]
	
	d.othr = rbind(d.othr, data.table(tag = tag,  
																		idx = rare,  
																		pos = POS[rare],
																		frq = AC[rare]))
	
	cat("\n")
}


splt = split(d.path, d.path$tag)
splt = lapply(splt, function(x) {
	state = as.numeric(factor(x$state))
	prev = 1
	seg0 = rep(T, nrow(x))
	seg1 = rep(T, nrow(x))
	for (curr in 2:length(state)) {
		flag = (state[prev] != state[curr])
		seg0[curr] = flag
		seg1[prev] = flag
		prev = curr
	}
	
	seq0 = which(seg0)
	seq1 = which(seg1)
	
	d = x[seg0]
	d$end.idx = x$idx[seg1]
	d$end.pos = x$pos[seg1]
	d
})
splt = rbindlist(splt)

d.viterbi.full.ibd = splt[which(splt$state == "IBD")]
d.viterbi.full.non = splt[which(splt$state == "NON")]


splt = split(d.halt, d.halt$tag)
splt = lapply(splt, function(x) {
	state = as.numeric(factor(x$state))
	prev = 1
	seg0 = rep(T, nrow(x))
	seg1 = rep(T, nrow(x))
	for (curr in 2:length(state)) {
		flag = (state[prev] != state[curr])
		seg0[curr] = flag
		seg1[prev] = flag
		prev = curr
	}
	
	seq0 = which(seg0)
	seq1 = which(seg1)
	
	d = x[seg0]
	d$end.idx = x$idx[seg1]
	d$end.pos = x$pos[seg1]
	d
})
splt = rbindlist(splt)

d.viterbi.halt.ibd = splt[which(splt$state == "IBD")]
d.viterbi.halt.non = splt[which(splt$state == "NON")]


d.outp = d.emis
d.outp$brk = cut(d.outp$pos, 5000)
d.outp = split(d.outp, list(d.outp$tag, d.outp$brk))
d.outp = lapply(d.outp, function(x) {
	z = x[1, ]
	#z$diff = mean(x$diff)
	z$diff = mean(x$non) - mean(x$ibd)
	#z$ibd = mean(x$ibd)
	#z$non = mean(x$non)
	z$idx.end = x$idx[nrow(x)]
	z$pos.end = x$pos[nrow(x)]
	z
})
d.outp = rbindlist(d.outp)

outp.res.lim = c(-0.055, 0.055)


gg = ggplot(d.post) + 
	facet_grid(tag~.) +
	geom_hline(yintercept = c(0, 0.5, 1), colour = "grey80") +
	geom_line(aes(x = pos, y = pprob)) +
	geom_rect(data = d.outp, aes(xmin = pos, xmax = pos.end, ymin = 1.1, ymax = 2.1, fill = diff)) + 
	#geom_raster(data = d.path, aes(x = pos, y = 1.6, fill = state)) +
	#geom_rect(data = d.segm, aes(xmin = pos, xmax = end.pos, ymin = 1.1, ymax = 2.1, fill = state)) + 
	geom_rect(data = d.viterbi.full.ibd, aes(xmin = pos, xmax = end.pos, ymin = 2.925, ymax = 3.6), fill = "orange") + 
	geom_rect(data = d.viterbi.full.non, aes(xmin = pos, xmax = end.pos, ymin = 2.925, ymax = 3.6), fill = "royalblue1") + 
	geom_rect(data = d.viterbi.halt.ibd, aes(xmin = pos, xmax = end.pos, ymin = 2.2, ymax = 2.875), fill = "orange") + 
	geom_rect(data = d.viterbi.halt.non, aes(xmin = pos, xmax = end.pos, ymin = 2.2, ymax = 2.875), fill = "royalblue1") + 
	geom_linerange(data = d.brks, aes(x = pos, ymin = 3.7, ymax = 4.2), colour = "grey60") +
	geom_linerange(data = d.true, aes(x = pos, ymin = 2.5, ymax = 3.3), colour="white") +
	geom_linerange(data = d.true, aes(x = pos, ymin = 2.5, ymax = 3.3), colour="black", linetype = "22") +
	geom_point(data = d.focl, aes(x = pos, y = 2.9), colour="black", size = 1.5) +
	geom_point(data = d.focl, aes(x = pos, y = 2.9), colour="white", size = 0.75) +
	geom_point(data = d.rare, aes(x = pos, y = 4.05), colour = "black", size = 1) +
	geom_point(data = d.othr, aes(x = pos, y = 3.85), colour = "black", size = 1, shape = 1) +
	#geom_point(data = d.post, aes(x = pos, y = 2.9)) +
	scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100e6, by = 2.5e6), labels = sprintf("%.1f", seq(0, 100e6, by = 2.5e6) / 1e6)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_fill_gradient2(low = "darkorange", high = "royalblue3", mid = "white", midpoint = 0, limits=outp.res.lim, expand=c(0,0), breaks=c(-0.05, -0.025, 0, 0.025, 0.05)) +
	#scale_fill_manual(values = c(IBD = "orange", NON = "royalblue1")) +
	scale_colour_manual(values = c(IBD = "orange", NON = "royalblue1")) +
	coord_cartesian(xlim = c(max(POS[1], fki$pos - 5e6), min(POS[length(POS)], fki$pos + 5e6))) +
	theme_few() +
	theme(aspect.ratio = 1/8,
		legend.title = element_blank(),
		legend.key.width = unit(0.05, "npc"),
		legend.position = "top") +
	xlab("Physical position (Mbp)") +
	ylab("Posterior probability") +
	ggtitle(sprintf("Focal rare variant: f%d, at position %.0f", AC[fki$index], fki$pos))


ggsave(gg, filename = sprintf("_plot.hmm.test_b.%d.%s.zoom.png", fki$index, prefix), width = 12, height = nrow(d.focl) * 1/2 + 3)





d.ibd = d.path[which(d.path$state == "IBD"), ]
d.non = d.path[which(d.path$state == "NON"), ]

ggplot() + 
	geom_point(data = d.ibd, aes(x = idx, y = 1), colour = "orange") + 
	geom_point(data = d.non, aes(x = idx, y = 0), colour = "royalblue1") +
	theme(aspect.ratio = 1/10)








