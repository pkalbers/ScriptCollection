#
# detect IBD with HMM
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]             # file prefix
prob.file = args[2]             # hmm probability file



### HMM methods


expected.age = function(k, n) {
	if (length(k) != 1) stop("No vector allowed")
	j = 2:n
	s = sum( choose(n - j, k - 1) * ((n - j + 1) / (n * (j - 1))) )
	2 * choose(n - 1, k)^-1 * s
}


get.init.prob = function(k, empirical.init) {
	p = empirical.init$ibd[ which(empirical.init$fk == k) ]
	c(IBD = p, NON = 1 - p)
}


hmm.viterbi = function (obs, init, trans, emiss) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts

	v = array(NA, c(n.sts, n.obs))

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





detect.shap.hmm = function(idx, g0, g1, G, POS, AC, hmm.init, hmm.emiss, sample.size, effect.size, recomb.rate) {

	foc = idx
	min = 1
	max = length(POS)

	g0 = G[, g0]
	g1 = G[, g1]

	if (g0[foc] != 1 || g0[foc] != 1) {
		return(c(NA, NA))
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
	lhs.init = get.init.prob(AC[foc], hmm.init)
	rhs.init = get.init.prob(AC[foc], hmm.init)


	# emiss probs
	lhs.emiss = hmm.emiss[ , , lhs.rng ]
	rhs.emiss = hmm.emiss[ , , rhs.rng ]


	# trans probs
	lhs.trans = array(c(NA, 0, NA, 1), c(2, 2, length(lhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL))
	rhs.trans = array(c(NA, 0, NA, 1), c(2, 2, length(rhs.obs)), dimnames = list(c("IBD", "NON"), c("IBD", "NON"), NULL))
	
	exp.age = expected.age(AC[foc], sample.size) * 2 * effect.size
	
	lhs.prb = exp( -2 * lhs.dist * recomb.rate * exp.age )
	rhs.prb = exp( -2 * rhs.dist * recomb.rate * exp.age )
	
	lhs.trans[1, 1, ] = lhs.prb
	lhs.trans[1, 2, ] = 1 - lhs.prb
	
	rhs.trans[1, 1, ] = rhs.prb
	rhs.trans[1, 2, ] = 1 - rhs.prb


	# HMM: Viterbi path
	lhs.path = hmm.viterbi(lhs.obs, lhs.init, lhs.trans, lhs.emiss)
	rhs.path = hmm.viterbi(rhs.obs, rhs.init, rhs.trans, rhs.emiss)


	lhs = which(rev(lhs.path) == 1)
	rhs = which(    rhs.path  == 1)

	if (length(lhs) == 0) { lhs = 0; } else { lhs = POS[ min(lhs) ]; }
	if (length(rhs) == 0) { rhs = 0; } else { rhs = POS[ max(rhs) + foc - 1 ]; }
	
	return(c(lhs, rhs))
}




### execute detections

run.hmm = function(pair, G, POS, AC, hmm.init, hmm.emiss.ibd, hmm.emiss.non, sample.size, effect.size, recomb.rate, log) {
	cat("###  HMM\n", file = log, append = T)

	# preparing emission probs
	hmm.emiss = array(NA, c(2, 6, length(POS)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))
	for (i in 1:length(POS)) {
		hmm.emiss[1 , , i] = hmm.emiss.ibd[AC[i] + 1, ]
		hmm.emiss[2 , , i] = hmm.emiss.non[AC[i] + 1, ]
	}
	
	pair$m.lhs = 0
	pair$m.rhs = 0

	for (i in 1:nrow(pair)) {
		if (i %% 100 == 0) cat(sprintf(" %d of %d\n", i, nrow(pair)), file = log, append = T)
		g.lr = detect.shap.hmm(pair$index[i], pair$g0[i], pair$g1[i], G, POS, AC, hmm.init, hmm.emiss, sample.size, effect.size, recomb.rate)
		pair$m.lhs[i] = g.lr[1]
		pair$m.rhs[i] = g.lr[2]
	}
	cat("\nOK\n", file = log, append = T)

	return(pair)
}





### run

load(sprintf("%s.RData", prefix))

load(prob.file)


G = sprintf("%s.G", prefix)
G = load.bigmatrix(G)
G = as.matrix(G)


AC = rowSums(G)

sample.size = 5000
effect.size = 10000
recomb.rate = 1e-8


setwd("pairs_hmm")


pair.files = dir(pattern = sprintf("^pairs_hmm\\.%s\\.[0-9]+\\.RData", prefix))



for (pair.file in sample(pair.files)) {

	cat("RUN: ", pair.file, "\n")

	run = sub("[^0-9]+([0-9]+)[^0-9]+", "\\1", pair.file)

	save.file = sprintf("result_hmm.%s",  pair.file)
	log.file  = sprintf("log_hmm.%s.txt", pair.file)

	Sys.sleep(runif(n = 1, min = 1/4, max = 2))

	if (file.exists(save.file) || file.exists(log.file)) {
		next
	}

	cat("", file = log.file, append = F)

	load(pair.file)

	if (file.exists(save.file)) next


	pair = run.hmm(pair, G, POS, AC, hmm.init, hmm.emiss.ibd, hmm.emiss.non, sample.size, effect.size, recomb.rate, log.file)


	if (file.exists(save.file)) next

	cat("Saving results\n", file = log.file, append = T)
	save(pair, file=save.file)
}






