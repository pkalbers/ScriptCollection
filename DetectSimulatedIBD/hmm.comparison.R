#
# HMM comparisons
#


### HMM functions:


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
		
		if (!is.null(cutoff) && z[i] < cutoff) {
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


### Other functions:

expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}

###


### HMM probability functions:

setup = list()
setup$init  = list()
setup$emiss = list()
setup$trans = list()


# Initial probs:

setup$init$Empirical = function(k, n, empirical.init) {
	p = empirical.init$ibd[ which(empirical.init$fk == k) ]
	c(IBD = p, NON = 1 - p)
}

setup$init$AgeDependent = function(k, n, empirical.init) {
	t = expected.age(k, n)
	p = (4/6) * (1 - exp(-6 * t))
	p = 1e-5
	c(IBD = p, NON = 1 - p)
}

setup$init$Fixed = function(k, n, empirical.init) {
	p = 1 - 1e-5
	c(IBD = p, NON = 1 - p)
}


# Emission probs:

setup$emiss$Empirical = function(k, empirical.ibd, empirical.non) {
	matrix(c(empirical.ibd[k + 1, ], empirical.non[k + 1, ]), nrow = 2, byrow = T, dimnames = list(c("IBD", "NON"), colnames(empirical.ibd)))
}

setup$emiss$Fixed = function(k, empirical.ibd, empirical.non) {
	matrix(c(
		0.7969959, 0.07396577, 0.0005191694, 0.05522923, 0.03675715, 0.03653279, 
		0.7703712, 0.10985430, 0.0181574900, 0.03687256, 0.03699387, 0.02775058
	), nrow = 2, byrow = T, dimnames = list(c("IBD", "NON"), c("00", "01", "02", "11", "12", "22")))
}


# Transision probs:

setup$trans$AgeDependent.left2right = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	t = expected.age(k, n)
	r = d * rec.rate
	
	q11 = exp(-4 * N * r * t)
	q10 = 1 - q11
	q01 = 0 
	q00 = 1 
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}

setup$trans$AgeDependent.ergodic = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	t = expected.age(k, n)
	r = d * rec.rate
	p = (4/6) * (1 - exp(-6 * t))
	
	q11 = exp(-4 * N * r * t)
	q10 = 1 - q11
	q01 = (p / (1 - p)) * q10
	q00 = 1 - q01
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}

setup$trans$Fixed.ergodic = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	p = 1e-5
	
	q10 = 5e-8
	q11 = 1 - q10
	q01 = (p / (1 - p)) * q10
	q00 = 1 - q01
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}

setup$trans$Fixed.left2right = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	q10 = 5e-8
	q11 = 1 - q10
	q01 = 0
	q00 = 1
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}




###


### Run:

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


prefix = "history.generror_1000g"
load(sprintf("%s.RData", prefix))
load("../truth/ibd_mrca.truth.RData")
load(sprintf("hmm.probs.%s.RData", prefix))

sample.size = 5000


# genotype data

G = as.matrix(load.bigmatrix(sprintf("%s.G", prefix)))

AC = rowSums(G)
AF = AC / (ncol(G) * 2)


# truth data

del = which(FKI$frq > 0.5)
if (length(del) > 0) {
	del = FKI$index[del]
	del = which(truth$index %in% del)
	if (length(del) > 0) {
		truth = truth[-del, ]
	}
}

del = which(! FKI$index %in% unique(truth$index))
if (length(del) > 0) {
	FKI = FKI[-del, ]
}

del = which(! truth$index %in% FKI$index)
if (length(del) > 0) {
	truth = truth[-del, ]
}

del = which(truth$wall.lhs | truth$wall.rhs)
if (length(del) > 0) {
	truth = truth[-del, ]
}


# common vars

dist.lhs = c(1, diff(POS))
dist.rhs = c(diff(POS), 1)

idx.min = 1
idx.max = length(POS)

pos.min = POS[1]
pos.max = POS[length(POS)]


splt = split(truth, truth$fk)
splt = lapply(splt, function(x, n) x[sample(nrow(x), size = min(n, nrow(x))), ], 10)
TRUTH = Reduce(rbind, splt)
TRUTH = TRUTH[sample(nrow(TRUTH)), ]

TRUTH$g0 = 0
i = which(TRUTH$h0 %% 2 != 0)
j = which(TRUTH$h0 %% 2 == 0)
TRUTH$g0[i] = (TRUTH$h0[i] + 1) / 2
TRUTH$g0[j] = TRUTH$h0[j] / 2

TRUTH$g1 = 0
i = which(TRUTH$h1 %% 2 != 0)
j = which(TRUTH$h1 %% 2 == 0)
TRUTH$g1[i] = (TRUTH$h1[i] + 1) / 2
TRUTH$g1[j] = TRUTH$h1[j] / 2





# tests

test = NULL

count = 0

for (select in 1:nrow(TRUTH)) {
	count = count + 1
	cat(sprintf(" %d of %d\n", count, nrow(TRUTH)))
	
	# select focal site
	idx.foc = TRUTH$index[select]
	idx.lhs = match(TRUTH$lhs[select], POS)
	idx.rhs = match(TRUTH$rhs[select], POS)
	
	pos.foc = POS[ TRUTH$index[select] ]
	pos.lhs = TRUTH$lhs[select]
	pos.rhs = TRUTH$rhs[select]
	
	selected.fk = TRUTH$fk[select]
	
	
	# genotypes
	selected.g0 = TRUTH$g0[select]
	selected.g1 = TRUTH$g1[select]
	
	g0 = G[, selected.g0]
	g1 = G[, selected.g1]
	
	obs.all = sprintf("%d%d", g0, g1)
	i = which(obs.all == "10"); if (length(i) > 0) obs.all[i] = "01"
	i = which(obs.all == "20"); if (length(i) > 0) obs.all[i] = "02"
	i = which(obs.all == "21"); if (length(i) > 0) obs.all[i] = "12"
	
	lhs.rng = idx.foc : idx.min
	rhs.rng = idx.foc : idx.max
	
	lhs.obs = obs.all[ lhs.rng ]
	rhs.obs = obs.all[ rhs.rng ]
	
	
	for (test.init in names(setup$init)) {
		
		# init probs
		init = setup$init[[ test.init ]](AC[idx.foc], sample.size, hmm.init)
		
		for (test.emiss in names(setup$emiss)) {
			
			# emiss probs
			lhs.emiss = sapply(lhs.rng, function(x) setup$emiss[[ test.emiss ]](AC[x], hmm.emiss.ibd, hmm.emiss.non), simplify = "array")
			rhs.emiss = sapply(rhs.rng, function(x) setup$emiss[[ test.emiss ]](AC[x], hmm.emiss.ibd, hmm.emiss.non), simplify = "array")
			
			for (test.trans in names(setup$trans)) {
				#cat(sprintf("\n### Init:  %s \n### Emiss: %s \n### Trans: %s \n", test.init, test.emiss, test.trans))
				
				# trans probs
				lhs.trans = sapply(lhs.rng, function(x) setup$trans[[ test.trans ]](dist.lhs[x], AC[idx.foc], sample.size), simplify = "array")
				rhs.trans = sapply(rhs.rng, function(x) setup$trans[[ test.trans ]](dist.rhs[x], AC[idx.foc], sample.size), simplify = "array")
				
				
				# HMM: Viterbi paths
				
				lhs.hmm.full = hmm.viterbi(lhs.obs, init, lhs.emiss, lhs.trans)
				rhs.hmm.full = hmm.viterbi(rhs.obs, init, rhs.emiss, rhs.trans)
				
				lhs.hmm.halt = hmm.viterbi(lhs.obs, init, lhs.emiss, lhs.trans, cutoff = 1e-16)
				rhs.hmm.halt = hmm.viterbi(rhs.obs, init, rhs.emiss, rhs.trans, cutoff = 1e-16)
				
				
				idx.lhs.hmm.full = idx.foc - which.max(lhs.hmm.full) + 1
				idx.rhs.hmm.full = idx.foc + which.max(rhs.hmm.full) - 1
				
				idx.lhs.hmm.halt = idx.foc - which.max(lhs.hmm.halt) + 1
				idx.rhs.hmm.halt = idx.foc + which.max(rhs.hmm.halt) - 1
				
				
				# DHG: breakpoints
				idx.lhs.dhg = idx.foc - which.max(as.numeric(lhs.obs == "02")) + 1
				idx.rhs.dhg = idx.foc + which.max(as.numeric(rhs.obs == "02")) - 1
				
				
				# store
				tmp = data.table( init  = test.init,
													emiss = test.emiss,
													trans = test.trans,
													fk = selected.fk,
													g0 = selected.g0,
													g1 = selected.g1,
													foc.idx = idx.foc,
													true.lhs.idx = idx.lhs,
													true.rhs.idx = idx.rhs,
													hmmf.lhs.idx = idx.lhs.hmm.full,
													hmmf.rhs.idx = idx.rhs.hmm.full,
													hmmh.lhs.idx = idx.lhs.hmm.halt,
													hmmh.rhs.idx = idx.rhs.hmm.halt,
													dhgt.lhs.idx = idx.lhs.dhg,
													dhgt.rhs.idx = idx.rhs.dhg,
													foc.pos = pos.foc,
													true.lhs.pos = pos.lhs,
													true.rhs.pos = pos.rhs,
													hmmf.lhs.pos = POS[ idx.lhs.hmm.full ],
													hmmf.rhs.pos = POS[ idx.rhs.hmm.full ],
													hmmh.lhs.pos = POS[ idx.lhs.hmm.halt ],
													hmmh.rhs.pos = POS[ idx.rhs.hmm.halt ],
													dhgt.lhs.pos = POS[ idx.lhs.dhg ],
													dhgt.rhs.pos = POS[ idx.rhs.dhg ] )
				
				test = rbind(test, tmp)
				
			}
		}
	}
	
}


save(test, file = sprintf("hmm.comparison.%s.RData", paste(sample(c(LETTERS, letters, as.character(0:9)), 16, replace = T), collapse = "")))


stop()



#####

tmp = NULL
for (file in dir(pattern = "^hmm\\.comparison\\..+\\.RData$")) {
	load(file)
	tmp = rbind(tmp, test)
}
tmp = rbindlist(lapply(split(tmp, tmp$foc.idx, tmp$g0, tmp$g1), function(x) if (nrow(x) == 24) return(x) else return(NULL) ))
test = tmp



a = abs(test$true.lhs.idx - test$hmmf.lhs.idx)
b = abs(test$true.rhs.idx - test$hmmf.rhs.idx)
d = data.table(foc = c(test$foc.idx, test$foc.idx),
							 trans = c(test$trans, test$trans),
							 emiss = c(test$emiss, test$emiss),
							 init = c(test$init, test$init),
							 dist = c(a, b))

ggplot(d) + 
	facet_wrap(~trans, nrow = 1) + 
	geom_boxplot(aes(emiss, dist+1, fill = init), outlier.shape=1, alpha = 0.5) + 
	scale_y_log10(breaks = c(1, 1*10^(1:8)), labels = sprintf("%d", c(1, 1*10^(1:8))), minor_breaks = as.vector(sapply(1:9, function(x) x*10^(0:10)))) + 
	#coord_flip() + 
	theme_linedraw() +
	theme(panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_line(size = 0.25, colour = "grey"),
				panel.grid.minor.y = element_line(size = 0.25, colour = "grey"))





gg = ggplot(test) + 
	facet_grid(init~emiss~trans) + 
	geom_abline(intercept = c(0, 0), slope = c(-1, 1), colour = "grey80") +
	geom_point(aes(x = log10(abs(true.rhs.pos - foc.pos)), y = log10(abs(hmmf.rhs.pos - foc.pos)),  colour = init), alpha = 0.25, size = 0.75, shape=1) +
	geom_point(aes(x = -log10(abs(true.lhs.pos - foc.pos)), y = log10(abs(foc.pos - hmmf.lhs.pos)), colour = init), alpha = 0.25, size = 0.75, shape=1) +
	scale_x_continuous(expand = c(0, 0)) + 
	scale_y_continuous(expand = c(0, 0)) + 
	geom_vline(xintercept = 0) +
	theme_few() +
	theme(aspect.ratio = 0.5) +
	coord_cartesian(xlim = c(-8.1, 8.1), ylim = c(-0.1, 8.1)) +
	xlab("log10 distance to TRUE breakpoint") +
	ylab("log10 distance to DETECTED breakpoint")

ggsave(filename = "_plot.hmm.comparison.pdf", plot = gg, width = 12, height = 8)


gg = ggplot(test) + 
	facet_grid(init~emiss~trans) + 
	geom_abline(intercept = c(0, 0), slope = c(-1, 1), colour = "grey80") +
	geom_point(aes(x =  log10(abs(true.rhs.pos - foc.pos)), y = ((hmmf.rhs.pos - hmmh.rhs.pos)), colour = init), alpha = 0.5, size = 0.75, shape=1) +
	geom_point(aes(x = -log10(abs(true.lhs.pos - foc.pos)), y = ((hmmf.lhs.pos - hmmh.lhs.pos)), colour = init), alpha = 0.5, size = 0.75, shape=1) +
	scale_x_continuous(expand = c(0, 0)) + 
	scale_y_continuous(expand = c(0, 0)) + 
	geom_vline(xintercept = 0) +
	theme_few() +
	theme(aspect.ratio = 0.5) +
	coord_cartesian(xlim = c(-8.1, 8.1)) + #, ylim = c(-0.1, 8.1)) +
	xlab("log10 distance, complete HMM") +
	ylab("Delta of complete - stopped HMM")

ggsave(filename = "_plot.hmm.comparison.hmmf-hmmh.pdf", plot = gg, width = 12, height = 8)





msr = function(tru, det, grp) {
	r = (tru - det)^2
	m = sapply(split(r, grp), mean)
	sqrt(m)
}

se <- function(x) sqrt(var(x)/length(x))


key = sprintf("I: %12s  |  E: %9s  |  T: %23s", test$init, test$emiss, test$trans)

z = msr(c(test$true.lhs.idx, test$true.rhs.idx, test$true.lhs.idx, test$true.rhs.idx, test$true.lhs.idx, test$true.rhs.idx),  
				c(test$hmmf.lhs.idx, test$hmmf.rhs.idx, test$hmmh.lhs.idx, test$hmmh.rhs.idx, test$dhgt.lhs.idx, test$dhgt.rhs.idx), 
				c(sprintf("%s_HMM, complete", key), sprintf("%s_HMM, complete", key), sprintf("%s_HMM, halted", key), sprintf("%s_HMM, halted", key), sprintf("%s_DHG", key), sprintf("%s_DHG", key)))

qq = qplot(x = z, y = sub("^(.+)_(.+)$", "\\1", names(z)), colour = sub("^(.+)_(.+)$", "\\2", names(z)),
					 xlab = "RMSE", ylab = "") + 
	theme(axis.text.y=element_text(family = "mono"), legend.title=element_blank())



z = msr(c(test$true.lhs.idx, test$true.rhs.idx),  
				c(test$hmmf.lhs.idx, test$hmmf.rhs.idx), 
				c(sprintf("%s_HMM, complete", key), sprintf("%s_HMM, complete", key)))

qq = qplot(x = z, y = sub("^(.+)_(.+)$", "\\1", names(z)), xlab = "RMSE", ylab = "") + 
	theme(axis.text.y=element_text(family = "mono"), legend.title=element_blank())


ggsave(filename = "_plot.hmm.comparison.RMSE.pdf", plot = qq, width = 12, height = 8)








