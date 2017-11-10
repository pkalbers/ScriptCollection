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
	if (length(k) != 1) stop("No vector allowed")
	j = 2:n
	s = sum( choose(n - j, k - 1) * ((n - j + 1) / (n * (j - 1))) )
	2 * choose(n - 1, k)^-1 * s
}

# expected.age = function(k, n) {
# 	x = k / n
# 	if (x == 1) return(2)
# 	((-2 * x) / (1 - x)) * log(x)
# }


get.trans.prob = function(d, k = NULL, n = NULL, t = NULL, N = 10000, rec.rate = 1e-8) {
	if (is.null(t)) {
		t = expected.age(k, n) * 2 * N
	}
	r = d * rec.rate
	p = exp(-2 * r * t)
	matrix(c(p, 0, 1 - p, 1), nrow = 2, ncol = 2, dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}


get.emiss.prob = function(k, empirical.ibd, empirical.non) {
	matrix(c(empirical.ibd[k + 1, ], empirical.non[k + 1, ]), nrow = 2, byrow = T, dimnames = list(c("IBD", "NON"), colnames(empirical.ibd)))
}


get.init.prob = function(k, empirical.init) {
	p = empirical.init$ibd[ which(empirical.init$fk == k) ]
	c(IBD = p, NON = 1 - p)
}



hmm.viterbi = function (obs, init, trans, emiss) {
	# states = names(init)
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	v = array(NA, c(n.sts, n.obs))
	#dimnames(v) = list(states, NULL)
	
	for (state in states) {
		v[state, 1] = log(init[state] * emiss[state, obs[1], 1])
	}
	
	for (i in 2:n.obs) {
		for (state in states) {
			max = NULL
			for (prev in states) {
				tmp = v[prev, i - 1] + log(trans[prev, state, i])
				max = max(max, tmp)
			}
			v[state, i] = log(emiss[state, obs[i], i]) + max
		}
	}
	
	path = rep(NA, n.obs)
	
	for (state in states) {
		if (max(v[, n.obs]) == v[state, n.obs]) {
			path[n.obs] = state
			break
		}
	}
	
	for (i in (n.obs - 1):1) {
		for (state in states) {
			p = path[i + 1]
			x = max(v[     , i] + log(trans[     , p, i]))
			y =     v[state, i] + log(trans[state, p, i])
			if (x == y) {
				path[i] = state
				break
			}
		}
	}
	
	return(path)
}



hmm.posterior = function (obs, init, trans, emiss) {
	# states = names(init)
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	forw = hmm.forward(obs, init, trans, emiss)
	back = hmm.backward(obs, init, trans, emiss)
	
	prob = forw[1, n.obs]
	
	for (i in 2:n.sts) {
		j = forw[i, n.obs]
		if (j > -Inf) {
			prob = j + log(1 + exp(prob - j))
		}
	}
	
	post = exp((forw + back) - prob)
	
	return(post)
}


hmm.forward = function (obs, init, trans, emiss) {
	# states = names(init)
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	forw = array(NA, c(n.sts, n.obs))
	#dimnames(forw) = list(states, NULL)
	
	for (state in states) {
		forw[state, 1] = log(init[state] * emiss[state, obs[1], 1])
	}
	
	for (i in 2:n.obs) {
		for (state in states) {
			lsum = -Inf
			for (prev.state in states) {
				tmp = forw[prev.state, i - 1] + log(trans[prev.state, state, i])
				if (tmp > -Inf) {
					lsum = tmp + log(1 + exp(lsum - tmp))
				}
			}
			forw[state, i] = log(emiss[state, obs[i], i]) + lsum
		}
	}
	
	return(forw)
}	


hmm.backward = function (obs, init, trans, emiss) {
	# states = names(init)
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	back = array(NA, c(n.sts, n.obs))
	#dimnames(back) = list(states, NULL)
	
	for (state in states) {
		back[state, n.obs] = log(1)
	}
	
	for (i in (n.obs - 1):1) {
		for (state in states) {
			lsum = -Inf
			for (next.state in states) {
				tmp = back[next.state, i + 1] + log(trans[state, next.state, i + 1] * emiss[next.state, obs[i + 1], i + 1])
				if (tmp > -Inf) {
					lsum = tmp + log(1 + exp(lsum - tmp))
				}
			}
			back[state, i] = lsum
		}
	}
	
	return(back)
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
chr.emiss = array(NA, c(2, 6, length(POS)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
for (i in 1:length(POS)) {
	chr.emiss[1 , , i] = hmm.emiss.ibd[AC[i] + 1, ]
	chr.emiss[2 , , i] = hmm.emiss.non[AC[i] + 1, ]
}



# random test

fki = FKI[sample(which(FKI$n.sharer == 8), 1), ]

g.sharers = as.numeric(strsplit(fki$g.sharer, "|", T)[[1]])
g.sharers = combn(g.sharers, 2)

h.sharers = as.numeric(strsplit(fki$h.sharer, "|", T)[[1]])
h.sharers = combn(h.sharers, 2)


d.brks = NULL
d.true = NULL
d.rare = NULL
d.path = NULL
d.post = NULL
d.focl = NULL

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
	
	
	# fetch truth
	
	tru = which(truth$index == foc & truth$h0 == pair.h[1] & truth$h1 == pair.h[2])
	
	if (length(tru) == 1) {
		true.position = c(truth$lhs[tru], truth$rhs[tru])
		true.index = match(true.position, POS)
		
		d.true = rbind(d.true, data.table(tag = paste(pair.g, collapse = " + "),  
																			idx = true.index,  
																			pos = true.position))
		
		cat("True LHS rel. index:", foc - true.index[1], "\n")
		cat("True RHS rel. index:", true.index[2] - foc, "\n")
		
		rare = truth$index[ which(truth$h0 == truth$h0[tru] & truth$h1 == truth$h1[tru]) ]
		
		d.rare = rbind(d.rare, data.table(tag = paste(pair.g, collapse = " + "),  
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
	
	lhs.dist = c(0, abs(diff(POS[ lhs.rng ])) )
	rhs.dist = c(0,     diff(POS[ rhs.rng ])  )
	
	lhs.freq = AC[ lhs.rng ]
	rhs.freq = AC[ rhs.rng ]
	
	
	# init probs
	cat("Getting initial probabilities ... \n")
	
	lhs.init = get.init.prob(AC[foc], hmm.init)
	rhs.init = get.init.prob(AC[foc], hmm.init)
	
	
	# emiss probs
	cat("Getting emission probabilities ... \n")
	
	lhs.emiss = chr.emiss[ , , lhs.rng ] # array(NA, c(2, 6, length(lhs.obs)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
	rhs.emiss = chr.emiss[ , , rhs.rng ] # array(NA, c(2, 6, length(rhs.obs)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
	
# 	for (i in 1:length(lhs.obs)) {
# 		lhs.emiss[ , , i] = get.emiss.prob(k = lhs.freq[i], hmm.emiss.ibd, hmm.emiss.non)
# 	}
# 	for (i in 1:length(rhs.obs)) {
# 		rhs.emiss[ , , i] = get.emiss.prob(k = rhs.freq[i], hmm.emiss.ibd, hmm.emiss.non)
# 	}
	
	
	# trans probs
	cat("Getting transition probabilities ... \n")
	
	lhs.trans = array(c(NA, 0, NA, 1), c(2, 2, length(lhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL)) # array(NA, c(2, 2, length(lhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL))
	rhs.trans = array(c(NA, 0, NA, 1), c(2, 2, length(rhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL)) # array(NA, c(2, 2, length(rhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL))
	
	exp.age = expected.age(AC[foc], sample.size) * 2 * effect.size
	
	lhs.prb = exp( -2 * lhs.dist * recomb.rate * exp.age )
	rhs.prb = exp( -2 * rhs.dist * recomb.rate * exp.age )
	
	lhs.trans[1, 1, ] = lhs.prb
	lhs.trans[1, 2, ] = 1 - lhs.prb
	
	rhs.trans[1, 1, ] = rhs.prb
	rhs.trans[1, 2, ] = 1 - rhs.prb
	
# 	for (i in 1:length(lhs.obs)) {
# 		lhs.trans[ , , i] = get.trans.prob(d = lhs.dist[i], t = exp.age)
# 	}
# 	for (i in 1:length(rhs.obs)) {
# 		rhs.trans[ , , i] = get.trans.prob(d = rhs.dist[i], t = exp.age)
# 	}
	
	
	# HMM: Viterbi path
	cat("Running Viterbi ... \n")
	lhs.path = hmm.viterbi(lhs.obs, lhs.init, lhs.trans, lhs.emiss)
	rhs.path = hmm.viterbi(rhs.obs, rhs.init, rhs.trans, rhs.emiss)
	
	
	# HMM: Posteriors
	cat("Running posteriors ... \n")
	
	lhs.trans[2, , ] = c(1e-8, 1 - 1e-8) ### required non-zero transitions
	rhs.trans[2, , ] = c(1e-8, 1 - 1e-8) ### required non-zero transitions
	
	lhs.post = hmm.posterior(lhs.obs, lhs.init, lhs.trans, lhs.emiss)
	rhs.post = hmm.posterior(rhs.obs, rhs.init, rhs.trans, rhs.emiss)
	
	
	path = rep("NON", length(POS))
	
	i = which(rev(lhs.path) == 1) # "IBD")
	if (length(i) > 0) {
		path[i] = "IBD"
	}
	i = which(rhs.path == 1) # "IBD")
	if (length(i) > 0) {
		path[i + foc - 1] = "IBD"
	}
	
	post = rep(NA, length(POS))
	post[ lhs.rng ] = lhs.post[1, ]
	post[ rhs.rng ] = rhs.post[1, ]
	
	inf = which(post == Inf)
	if (length(inf) > 0) {
		post[inf] = 1
	}
	
	brks = which(obs.all == "02")
	
	tag = paste(pair.g, collapse = " + ")
	
	d.brks = rbind(d.brks, data.table(tag = tag,  idx = brks,           pos = POS[brks]))
	d.path = rbind(d.path, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,         state = path))
	d.post = rbind(d.post, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,         pprob = post))
	d.focl = rbind(d.focl, data.table(tag = tag,  idx = foc,            pos = POS[foc]))
	
	cat("\n")
}




gg = ggplot(d.path) + 
	facet_grid(tag~.) +
	geom_raster(aes(x = idx, y = 1.6, fill = state)) +
	geom_hline(yintercept = c(0, 0.5, 1), colour = "grey80") +
	geom_line(data = d.post, aes(x = idx, y = pprob)) +
	geom_linerange(data = d.brks, aes(x = idx, ymin = 2.2, ymax = 2.7), colour = "grey60") +
	geom_linerange(data = d.true, aes(x = idx, ymin = 1.1, ymax = 2.1), colour="white") +
	geom_linerange(data = d.true, aes(x = idx, ymin = 1.1, ymax = 2.1), colour="black", linetype = "22") +
	geom_point(data = d.focl, aes(x = idx, y = 1.6), colour="black", size = 1.5) +
	geom_point(data = d.focl, aes(x = idx, y = 1.6), colour="white", size = 0.75) +
	geom_point(data = d.rare, aes(x = idx, y = 2.45), colour = "black", size = 1) +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_fill_manual(values = c(IBD = "orange", NON = "royalblue1")) +
	theme_few() +
	theme(aspect.ratio = 1/10,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Chromosome position (index)") +
	ylab("Posterior probability")

ggsave(gg, filename = sprintf("_plot.hmm.test.%d.%s.pdf", fki$index, prefix), width = 12, height = nrow(d.true) * 1/2 + 3)






