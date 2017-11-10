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

# expected.age = function(k, n) {
# 	x = k / n
# 	if (x == 0) return(0)
# 	if (x == 1) return(2)
# 	((-2 * x) / (1 - x)) * log(x)
# }


get.init.prob = function(k, n) {
	t = expected.age(k, n)
	p = (3/6) * (1 - exp(-6 * t))
	p = 1e-5
	c(IBD = p, NON = 1 - p)
}


get.trans.prob = function(d, k = NULL, n = NULL, N = 10000, rec.rate = 1e-8) {
	t = expected.age(k, n)
	r = d * rec.rate
	p = (3/6) * (1 - exp(-6 * t))
	
	q11 = exp(-4 * N * r * t)
	q10 = 1 - q11
	q01 = (p / (1 - p)) * q10
	q00 = 1 - q01
	
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
for (d in c(0, 5, 10, 50, 100, 500, 1000, seq(10000, 1e6, by = 1000))) {
	for (k in c(2, 5, 10, 15, 20, 25, 5000-1)) {
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


hmm.viterbi = function (obs, init, trans, emiss) {
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
		for (state in states) {
			max = 0
			for (prev.state in states) {
				tmp = v[prev.state, i - 1] * trans[prev.state, state, i]
				max = max(max, tmp)
			}
			v[state, i] = emiss[state, obs[i], i] * max
		}
		
		w[i] = max(v[, i])
		v[, i] = v[, i] / w[i]
	}
	
	path = rep(NA, n.obs)
	
	tmp = which.max(v[, n.obs])
	path[i] = states[tmp]
	# 	for (state in states) {
	# 		if (max(v[, n.obs]) == v[state, n.obs]) {
	# 			path[n.obs] = state
	# 			break
	# 		}
	# 	}
	
	for (i in (n.obs - 1):1) {
		tmp = which.max(v[ , i] * trans[ , path[i + 1], i])
		path[i] = states[tmp]
		# 		max = max(v[ , i] * trans[ , path[i + 1], i])
		# 		for (state in states) {
		# 			if (max == v[state, i] * trans[state, path[i + 1], i]) {
		# 				path[i] = state
		# 				break
		# 			}
		# 		}
	}
	
	return(list(path = path, w = w))
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


hmm.loglikelihood = function (obs, init, emiss, trans) {
	n.obs = length(obs)
	n.sts = length(init)
	states = 1:n.sts
	
	fwd = hmm.fwd(obs, init, emiss, trans)
	bwd = hmm.bwd(obs, init, emiss, trans)
	
	logl.state = array(NA, c(n.sts, n.obs))
	logl.model = array(0, n.obs)
	
	wght.fwp = cumsum(log(fwd$w))
	wght.bwp = rev(cumsum(rev(log(bwd$w))))
	
	for (state in states) {
		logl.state[state, ] = log(fwd$f[state, ] * bwd$b[state, ]) + wght.fwp + wght.bwp
		logl.model = logl.model + (fwd$f[state, ] * bwd$b[state, ])
	}
	
	logl.model = log(logl.model) + wght.fwp + wght.bwp
	
	logp.fwd = wght.fwp[n.obs] + log(fwd$f[1, n.obs])
	logp.bwd = wght.bwp[  1  ] + log(bwd$b[1,   1  ] * emiss[1, obs[  1  ],   1  ] * init[1])
	
	for (i in 2:n.sts) {
		logp.fwd = logp.fwd + log(1 + exp((wght.fwp[n.obs] + log(fwd$f[i, n.obs])) - logp.fwd))
		logp.bwd = logp.bwd + log(1 + exp((wght.bwp[  1  ] + log(bwd$b[i,   1  ] * emiss[i, obs[  1  ],   1  ] * init[i])) - logp.bwd))
	}
	
	return(list(logl.state = logl.state, 
							logl.model = logl.model,
							logp.fwd = as.numeric(logp.fwd), 
							logp.bwd = as.numeric(logp.bwd)))
}




###


prefix = "1000GP_Phase3_chr20"


load("../generror_1000g/hmm.probs.history.generror_1000g.RData")

load(sprintf("%s.RData", prefix))
load(sprintf("%s.G.RData", prefix))


panel = read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)


AC = rowSums(G)
AF = AC / (ncol(G) * 2)



cross = sort(unique(AF)) * (nrow(hmm.emiss.ibd) - 1)
count = 0:(nrow(hmm.emiss.ibd) - 1)

approx.hmm.emiss.ibd = array(NA, c(length(cross), 6), dimnames = list(NULL, colnames(hmm.emiss.ibd)))
approx.hmm.emiss.non = array(NA, c(length(cross), 6), dimnames = list(NULL, colnames(hmm.emiss.non)))

for (pair in colnames(hmm.emiss.ibd)) {
	approx.hmm.emiss.ibd[, pair] = approx(count, hmm.emiss.ibd[, pair], cross)$y
	approx.hmm.emiss.non[, pair] = approx(count, hmm.emiss.non[, pair], cross)$y
}



emiss = array(NA, c(2, 6, length(POS)), dimnames = list(c("IBD", "NON"), colnames(hmm.emiss.ibd), NULL))

for (i in 1:length(POS)) {
	tmp1 = approx.hmm.emiss.ibd[AC[i] + 1, ] + 1e-8
	tmp2 = approx.hmm.emiss.non[AC[i] + 1, ] + 1e-8
	
	emiss[1 , , i] = tmp1 / sum(tmp1)
	emiss[2 , , i] = tmp2 / sum(tmp2)
}


# random test

rare.vars = which(AC >=2 & AC <= 25)

fki = list()
fki$index = sample(which(AC == 4), 1)
fki$pos = POS[fki$index]
fki$frq = AF[fki$index]



g.sharers = which(G[fki$index, ] == 1)
g.sharers = combn(g.sharers, 2)


init = get.init.prob(AC[fki$index], (ncol(G) * 2))


dist.fwd = c(diff(POS), 1)
dist.bwd = c(1, diff(POS))

trans = sapply(1:length(POS), function(x) get.trans.prob(dist.fwd[x], AC[fki$index], (ncol(G) * 2)), simplify = "array")



d.brks = NULL
d.rare = NULL
d.path = NULL
d.post = NULL
d.focl = NULL

for (share in 1:ncol(g.sharers)) {
	
	pair.g = g.sharers[, share]
	
	g0 = G[, pair.g[1]]
	g1 = G[, pair.g[2]]
	
	foc = fki$index
	
	cat("fk =", AC[foc], "\n")
	
	if (g0[foc] != 1 || g0[foc] != 1) {
		cat("Rare variant not detected!\n")
		next
	}
	
	
	# observed genotye pairs
	
	obs = sprintf("%d%d", g0, g1)
	i = which(obs == "10"); if (length(i) > 0) obs[i] = "01"
	i = which(obs == "20"); if (length(i) > 0) obs[i] = "02"
	i = which(obs == "21"); if (length(i) > 0) obs[i] = "12"
	
	
	# HMM: Viterbi path
	cat("Running Viterbi ... \n")
	path = hmm.viterbi(obs, init, trans, emiss)
	viterbi.path = path$path
	path.weights = path$w
	
	i = which(viterbi.path == 1); if (length(i) > 0) viterbi.path[i] = "IBD"
	i = which(viterbi.path == 2); if (length(i) > 0) viterbi.path[i] = "NON"
	
	
	# HMM: Posteriors
	cat("Running posteriors ... \n")
	logl = hmm.loglikelihood(obs, init, emiss, trans)
	post = exp(logl$logl.state - logl$logp.fwd)
	#post = exp(logl$logl.state - (logl$logl.state[1, ] + log(1 + exp(logl$logl.state[2, ] - logl$logl.state[1, ]))))
	
	
	brks = which(obs == "02")
	
	tag = paste(panel$sample[pair.g], collapse = " + ")
	
	d.brks = rbind(d.brks, data.table(tag = tag,  idx = brks,           pos = POS[brks]))
	d.path = rbind(d.path, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,         state = viterbi.path))
	d.post = rbind(d.post, data.table(tag = tag,  idx = 1:length(POS),  pos = POS,         pprob = post[1, ]))
	d.focl = rbind(d.focl, data.table(tag = tag,  idx = foc,            pos = POS[foc]))
	
	# other rare vars shared
	is.shared = which(g0[rare.vars] == 1 & g1[rare.vars] == 1)
	rare = rare.vars[is.shared]
	
	d.rare = rbind(d.rare, data.table(tag = tag,  
																		idx = rare,  
																		pos = POS[rare],
																		frq = AC[rare]))
	
	cat("\n")
}


d.segm = split(d.path, d.path$tag)
d.segm = lapply(d.segm, function(x) {
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
d.segm = rbindlist(d.segm)


gg = ggplot(d.post) + 
	facet_grid(tag~.) +
	geom_hline(yintercept = c(0, 0.5, 1), colour = "grey80") +
	geom_line(aes(x = pos, y = pprob)) +
	#geom_raster(data = d.path, aes(x = pos, y = 1.6, fill = state)) +
	geom_rect(data = d.segm, aes(xmin = pos, xmax = end.pos, ymin = 1.1, ymax = 2.1, fill = state)) + 
	geom_linerange(data = d.brks, aes(x = pos, ymin = 2.2, ymax = 2.7), colour = "grey60") +
	geom_point(data = d.focl, aes(x = pos, y = 1.6), colour="black", size = 1.5) +
	geom_point(data = d.focl, aes(x = pos, y = 1.6), colour="white", size = 0.75) +
	geom_point(data = d.rare, aes(x = pos, y = 2.45), colour = "black", size = 1, shape = 1) +
	#geom_point(data = d.post, aes(x = pos, y = 2.9)) +
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_continuous(breaks = c(0, 0.5, 1)) +
	scale_fill_manual(values = c(IBD = "orange", NON = "royalblue1")) +
	coord_cartesian(xlim = fki$pos + c(-5e6, 5e6)) +
	theme_few() +
	theme(aspect.ratio = 1/15,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Physical position") +
	ylab("Posterior probability") +
	ggtitle(sprintf("Focal rare variant: f%d, at position %d", AC[fki$index], fki$pos))
#gg
ggsave(gg, filename = sprintf("_plot.hmm.test_b.%d.%s.zoom.png", fki$index, prefix), width = 12, height = nrow(d.focl) * 1/2 + 3)






