#
# HMM to detect shared haplotype on genotypes
#

library(data.table)

args = commandArgs(T)

prefix = args[1] # file prefix


### HMM functions

make.genotype.pair = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}


expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}


get.prob.trans = function(r, k, n , N) {
	t = expected.age(k, n)
	
	q11 = exp(-4 * N * r * t)
	q10 = 1 - q11
	q01 = 0 
	q00 = 1 
	
	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}


hmm.viterbi = function (obs, init, emiss, trans) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	v = array(NA, c(n.sts, n.obs))
	w = array(NA, n.obs)
	
	for (state in states) {
		v[state, 1] = init[state] * emiss[state, obs[1], 1]
	}
	
	w[1] = max(v[, 1])
	v[, 1] = v[, 1] / w[1]
	
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
		
		w[i] = max(v[, i])
		v[, i] = v[, i] / w[i]
	}
	
	path = rep(NA, n.obs)
	
	tmp = which.max(v[, n.obs])
	path[i] = states[tmp]
	
	for (i in (n.obs - 1):1) {
		tmp = which.max(v[ , i] * trans[ , path[i + 1], i])
		path[i] = states[tmp]
	}
	
	return(path)
}




### detection methods

detect.shap.hmm = function(focal, G, POS, INIT, EMISS, TRANS) {
	foc = focal$index
	fk  = focal$fk
	g0  = focal$g0
	g1  = focal$g1
	
	fkc = as.character(fk)
	
	lhs.rng = foc : 1
	rhs.rng = foc : (length(POS))
	
	# observations
	obs = make.genotype.pair(G[, g0], G[, g1])
	lhs.obs = obs[ lhs.rng ]
	rhs.obs = obs[ rhs.rng ]
	
	# initial probs
	init = INIT[, fkc]
	
	# emission probs
	lhs.emiss = EMISS[ , , lhs.rng]
	rhs.emiss = EMISS[ , , rhs.rng]
	
	# transition probs
	lhs.trans = TRANS$LHS[[ fkc ]][ , , lhs.rng]
	rhs.trans = TRANS$RHS[[ fkc ]][ , , rhs.rng]
	
	lhs.trans[ , , 1] = c(1, 0, 0, 1)
	rhs.trans[ , , 1] = c(1, 0, 0, 1)
	
	# Viterbi path decoding
	lhs.path = hmm.viterbi(lhs.obs, init, lhs.emiss, lhs.trans)
	rhs.path = hmm.viterbi(rhs.obs, init, rhs.emiss, rhs.trans)
	
	lhs = match(2, lhs.path)
	rhs = match(2, rhs.path); if (is.na(rhs)) rhs = length(rhs.path)
	
	if (is.na(lhs)) {
		lhs = POS[1]
	} else {
		if (lhs == 1) {
			lhs = POS[foc]
		} else {
			lhs = POS[ lhs.rng[lhs] + 1 ]
		}
	}
	
	if (is.na(rhs)) {
		rhs = POS[length(POS)]
	} else {
		if (rhs == 1) {
			rhs = POS[foc]
		} else {
			rhs = POS[ rhs.rng[rhs] - 1 ]
		}
	}
	
	return(c(lhs, rhs))
}


### execute detections

run.hmm = function(pair, G, POS, RATE, Ne, INIT, EMISS, TRANS, log.file, save.file) {
	cat("### HMM\n", file = log.file, append = T)
	
	pair$hmm.lhs = 0
	pair$hmm.rhs = 0
	
	for (i in 1:nrow(pair)) {
		hmm.lr = detect.shap.hmm(pair[i, ], G, POS, INIT, EMISS, TRANS)
		pair$hmm.lhs[i] = hmm.lr[1]
		pair$hmm.rhs[i] = hmm.lr[2]
		
		if (i %% 100 == 0) {
			cat(sprintf(" %d of %d\n", i, nrow(pair)), file = log.file, append = T)
			save(pair, file=save.file)
		}
	}
	cat("\nOK\n", file = log.file, append = T)
	
	return(pair)
}




### run

load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.G.RData", prefix)) # genotypes
load(sprintf("result.hmm_probs.%s.RData", prefix)) # HMM probs


Ne = 7300


lhs.dist = c(diff(POS), 0) * RATE * 1e-8
rhs.dist = c(0, diff(POS)) * RATE * 1e-8


setwd("pairs_hmm")


pair.files = dir(pattern = sprintf("^pairs\\.%s\\.[0-9]+\\.RData", prefix))

for (pair.file in sample(pair.files)) {
	
	cat("RUN: ", pair.file, "\n")
	
	save.file = sprintf("result_hmm.%s", pair.file)
	log.file  = sprintf("log_hmm.%s.txt", pair.file)
	tmp.file  = sprintf("tmp_hmm.%s", pair.file)
	
	Sys.sleep(runif(n = 1, min = 1/4, max = 2))
	
	if (file.exists(save.file) || file.exists(log.file)) {
		next
	}
	
	cat("", file = log.file, append = F)
	
	load(pair.file)
	
	if (file.exists(save.file)) next
	
	cat("Preparing transition matrices ...")
	trans = list(LHS = list(), RHS = list())
	for (fk in sort(unique(pair$fk))) {
		cat(sprintf(" %d ...", fk))
		fkc = as.character(fk)
		trans$LHS[[fkc]] = sapply(lhs.dist, get.prob.trans, fk, ncol(G) * 2, Ne, simplify = "array")
		trans$RHS[[fkc]] = sapply(rhs.dist, get.prob.trans, fk, ncol(G) * 2, Ne, simplify = "array")
	}
	cat(" OK\n")
	
	pair = run.hmm(pair, G, POS, RATE, Ne, emp.init, emp.emiss, trans, log.file, tmp.file)
	
	if (file.exists(save.file)) next
	
	cat("Saving results\n", file = log.file, append = T)
	save(pair, file=save.file)
}













