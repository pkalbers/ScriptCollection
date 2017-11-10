#
# HMM to detect shared haplotype on genotypes
#

library(data.table)
library(parallel)


### HMM functions

make.genotype.pair = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}

make.haplotype.pair = function(h0, h1, mask = c("00"="00", "01"="01", "10"="01", "11"="11")) {
	h = sprintf("%d%d", h0, h1)
	as.vector(mask[h])
}


expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}


get.prob.trans = function(r, k, n , N) {
	t = expected.age(k, n)

	q11 = exp(-2 * N * r * t)
	q10 = 1 - q11
	q01 = 0
	q00 = 1

	array(c(q11, q01, q10, q00), c(2, 2), dimnames = list(c("IBD", "NON"), c("IBD", "NON")))
}


get.prob.inits = function(k) {
	ibd = rep(1 - 1e-08, length(k))
	non = rep(1e-08, length(k))

	array(matrix(c(ibd, non), nrow = 2, byrow = T), c(2, length(k)), dimnames = list(c("IBD", "NON"), as.character(k)))
}


get.prob.emiss = function(q, p = 1 - q) {
	n00 = p^2
	n01 = 2*p*q
	n11 = q^2

	i00 = p
	i01 = 2*p*q * 0.0001  #
	i11 = q
	s = i00 + i01 + i11
	i00 = i00 / s
	i01 = i01 / s
	i11 = i11 / s

	array(matrix(c(i00, n00, i01, n01, i11, n11), nrow = 6, byrow = T), c(2, 3, length(q)), dimnames = list(c("IBD", "NON"), c("00", "01", "11"), NULL))
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

detect.shap.hmm = function(focal, INITS, EMISS, TRANS, n) {
	foc = focal$foc
	fkc = as.character(focal$fk)

	lhs.rng = foc : 1
	rhs.rng = foc : n

	# observations
	obs = make.haplotype.pair(focal$h0, focal$h1)
	lhs.obs = obs[ lhs.rng ]
	rhs.obs = obs[ rhs.rng ]

	# initial probs
	init = INITS[, fkc]

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

	lhs.decode = function(path, rng, foc) {
		k = match(2, path)
		if (is.na(k)) {
			k = 1
		} else {
			if (k == 1) {
				k = foc
				#if (foc != 1) k = foc - 1
			} else {
				k = rng[k] + 1
			}
		}
		k
	}
	
	rhs.decode = function(path, rng, foc, n) {
		k = match(2, path)
		if (is.na(k)) {
			k = n
		} else {
			if (k == 1) {
				k = foc
				#if (foc != n) k = foc + 1
			} else {
				k = rng[k]  - 1
			}
		}
		k
	}
	
	lhs = lhs.decode(lhs.path, lhs.rng, foc)
	rhs = rhs.decode(rhs.path, rhs.rng, foc, n)
	
	
	# correct phase
	LHS = lhs
	RHS = rhs
	
	while (LHS != 1) {
		h0 = focal$h0
		h1 = focal$h1
		x0 = focal$x0
		x1 = focal$x1
		
		rng = 1:LHS
		
		tmp = h0[rng]
		h0[rng] = x0[rng]
		x0[rng] = tmp
		
		tmp = h1[rng]
		h1[rng] = x1[rng]
		x1[rng] = tmp
		
		l0.obs = make.haplotype.pair(h0[ lhs.rng ], focal$h1[ lhs.rng ])
		l0.path = hmm.viterbi(l0.obs, init, lhs.emiss, lhs.trans)
		l0 = lhs.decode(l0.path, lhs.rng, foc)
		
		l1.obs = make.haplotype.pair(focal$h0[ lhs.rng ], h1[ lhs.rng ])
		l1.path = hmm.viterbi(l1.obs, init, lhs.emiss, lhs.trans)
		l1 = lhs.decode(l1.path, lhs.rng, foc)
		
		l2.obs = make.haplotype.pair(h0[ lhs.rng ], h1[ lhs.rng ])
		l2.path = hmm.viterbi(l2.obs, init, lhs.emiss, lhs.trans)
		l2 = lhs.decode(l2.path, lhs.rng, foc)
		
		l = c(l0, l1, l2)  ;  print(c(LHS, l))
		x = which.min(l)
		
		if (l[x] < LHS) {
			LHS = l[x]
			if (x != 2) { focal$h0 = h0; focal$x0 = x0; }
			if (x != 1) { focal$h1 = h1; focal$x1 = x1; }
		} else {
			break
		}
	}
	
	while (RHS != n) {
		h0 = focal$h0
		h1 = focal$h1
		x0 = focal$x0
		x1 = focal$x1
		
		rng = RHS:n
		
		tmp = h0[rng]
		h0[rng] = x0[rng]
		x0[rng] = tmp
		
		tmp = h1[rng]
		h1[rng] = x1[rng]
		x1[rng] = tmp
		
		r0.obs = make.haplotype.pair(h0[ rhs.rng ], focal$h1[ rhs.rng ])
		r0.path = hmm.viterbi(r0.obs, init, rhs.emiss, rhs.trans)
		r0 = rhs.decode(r0.path, rhs.rng, foc, n)
		
		r1.obs = make.haplotype.pair(focal$h0[ rhs.rng ], h1[ rhs.rng ])
		r1.path = hmm.viterbi(r1.obs, init, rhs.emiss, rhs.trans)
		r1 = rhs.decode(r1.path, rhs.rng, foc, n)
		
		r2.obs = make.haplotype.pair(h0[ rhs.rng ], h1[ rhs.rng ])
		r2.path = hmm.viterbi(r2.obs, init, rhs.emiss, rhs.trans)
		r2 = rhs.decode(r2.path, rhs.rng, foc, n)
		
		r = c(r0, r1, r2)  ;   print(c(RHS, r))
		x = which.max(r)
		
		if (r[x] > RHS) {
			RHS = r[x]
			if (x != 2) { focal$h0 = h0; focal$x0 = x0; }
			if (x != 1) { focal$h1 = h1; focal$x1 = x1; }
		} else {
			break
		}
	}

	return(c(LHS, RHS))
}


sub.detect.shap.hmm = function(focal, INITS, EMISS, TRANS, n) {
	mclapply(focal, detect.shap.hmm, INITS, EMISS, TRANS, n, mc.cores = length(focal))
}


run.detect.shap.hmm = function(focal, H, POS, INITS, EMISS, TRANS) {
	n = nrow(focal)

	L = list()
	for (i in 1:n) {
		x0 = focal$h0[i] + 1
		x1 = focal$h1[i] + 1
		if (focal$h0[i] %% 2 == 0) x0 = focal$h0[i] - 1
		if (focal$h1[i] %% 2 == 0) x1 = focal$h1[i] - 1
		L[[i]] = list(foc = focal$index[i],
									fk  = focal$fk[i],
									h0  = H[, focal$h0[i]],
									h1  = H[, focal$h1[i]],
									x0  = H[, x0],
									x1  = H[, x1])
	}

	idx = sub.detect.shap.hmm(L, INITS, EMISS, TRANS, nrow(H))

	lhs = sapply(idx, function(x) x[1])
	rhs = sapply(idx, function(x) x[2])


	####

	focal$cc.lhs.idx = lhs
	focal$cc.rhs.idx = rhs

	focal$cc.lhs.pos = POS[lhs]
	focal$cc.rhs.pos = POS[rhs]

	focal
}



######################


ncores = 100

#prefix = "OutOfAfricaHapMap20"
prefix = "OutOfAfricaHapMap20.GenErr_1000G"


load(sprintf("data.%s.RData", prefix))
#load(sprintf("data.%s.H.RData", prefix))
load(sprintf("data.%s.P.RData", prefix))
#load(sprintf("result.match_hmm.%s.RData", prefix))
load(sprintf("result.match_hmm_hap.%s.RData", prefix))


#HP = H
HP = P


Ne = 7300


lhs.dist = c(diff(POS), 0) * RATE * 1e-8
rhs.dist = c(0, diff(POS)) * RATE * 1e-8

fk.rng = c(2, 5, 10, 15, 20, 25) # 2:25

cat("Initialising ... ")

INITS = get.prob.inits(fk.rng)
EMISS = get.prob.emiss(AAF)
TRANS = list(LHS = list(), RHS = list())

for (fk in fk.rng) {
	cat(fk, " ")
	fkc = as.character(fk)
	TRANS$LHS[[fkc]] = sapply(lhs.dist, get.prob.trans, fk, ncol(HP), Ne, simplify = "array")  ###
	TRANS$RHS[[fkc]] = sapply(rhs.dist, get.prob.trans, fk, ncol(HP), Ne, simplify = "array")  ###
}

cat("OK \n")


# cat("Sampling ... ")
# 
# match = match[sample(1:nrow(match)), ]
# match = match[order(match$fk), ]
# 
# x = which(match$true.wall)
# if (length(x) > 0) match = match[-x, ]
# 
# x = sprintf("%d %d %d %d", match$h0, match$h1, match$true.lhs.idx, match$true.rhs.idx)
# x = which(duplicated(x))
# if (length(x) > 0) match = match[-x, ]
# 
# 
# x = sprintf("%d %d %.4f %.4f", match$h0, match$h1, match$h.lhs, match$h.rhs)
# x = which(duplicated(x))
# if (length(x) > 0) match = match[-x, ]
# 
# x = sprintf("%d %d %.4f %.4f", match$h0, match$h1, match$p.lhs, match$p.rhs)
# x = which(duplicated(x))
# if (length(x) > 0) match = match[-x, ]
# 
# x = sprintf("%d %d %.4f %.4f", match$h0, match$h1, match$g.lhs, match$g.rhs)
# x = which(duplicated(x))
# if (length(x) > 0) match = match[-x, ]
# 
# 
# x = which(match$fk %in% fk.rng)
# match = match[x, ]
# 
# 
# match = split(match, match$fk)
# match = lapply(match, function(x) {
# 	n = min(nrow(x), 500)
# 	i = sort(sample(1:nrow(x), n))
# 	x[i, ]
# })
# match = rbindlist(match)
# 
# save(match, file = sprintf("result.match_hmm_hap.%s.RData", prefix))
# 
# match = split(match, 1:nrow(match) %% ceiling(nrow(match)/ncores))
# 
# cat("OK \n")


match = split(match, 1:nrow(match) %% ceiling(nrow(match)/ncores))



cat("Running:\n")
for (tag in names(match)) {
	cat(tag, "of", length(match), "\n")
	focal = match[[tag]]
	match[[tag]] = run.detect.shap.hmm(focal, HP, POS, INITS, EMISS, TRANS)  ######
}
cat("DONE \n")


match = rbindlist(match)


save(match, file = sprintf("result.match_hmm_hap.%s.RData", prefix))








stop()

####################


library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20"
#prefix = "OutOfAfricaHapMap20.GenErr_1000G"


load(sprintf("data.%s.RData", prefix))
load(sprintf("result.match_hmm_hap.%s.RData", prefix))


d = rbind(data.table(side="LHS", fk = match$fk, type = "(c) DGT, genotypes",            det = match$pos - match$g.lhs + 1,   tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(a) FGT, true haplotypes ",     det = match$pos - match$h.lhs + 1,   tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(b) FGT, phased haplotypes ",   det = match$pos - match$p.lhs + 1,   tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(f) G-HMM, genotypes",          det = match$pos - match$hmm.lhs,     tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(d) H-HMM, true haplotypes",    det = match$pos - match$hh.lhs.pos,  tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(e) H-HMM, phased haplotypes",  det = match$pos - match$pp.lhs.pos,  tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(g) H-HMM*, true haplotypes",   det = match$pos - match$ch.lhs.pos,  tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="LHS", fk = match$fk, type = "(h) H-HMM*, phased haplotypes", det = match$pos - match$cc.lhs.pos,  tru = match$pos - POS[match$true.lhs.idx] + 1),
					data.table(side="RHS", fk = match$fk, type = "(c) DGT, genotypes",            det = match$g.rhs - match$pos + 1,   tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(a) FGT, true haplotypes ",     det = match$h.rhs - match$pos + 1,   tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(b) FGT, phased haplotypes ",   det = match$p.rhs - match$pos + 1,   tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(f) G-HMM, genotypes",          det = match$hmm.rhs - match$pos,     tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(d) H-HMM, true haplotypes",    det = match$hh.rhs.pos - match$pos,  tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(e) H-HMM, phased haplotypes",  det = match$pp.rhs.pos - match$pos,  tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(g) H-HMM*, true haplotypes",   det = match$ch.rhs.pos - match$pos,  tru = POS[match$true.rhs.idx] - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(h) H-HMM*, phased haplotypes", det = match$cc.rhs.pos - match$pos,  tru = POS[match$true.rhs.idx] - match$pos + 1)
)



d$map = (d$det / d$tru)

p = d


p = split(p, list(p$type, p$fk))
p = lapply(p, function(z) {
	x = seq(0, 3, by=0.01)
	y = ecdf(z$map)(x)
	data.table(type = z$type[1], fk = z$fk[1], x=x, y=y)
})
p = rbindlist(p)

p$fk = factor(p$fk, levels = c(2:25), ordered = T)

hist = ggplot(data=p) +
	facet_wrap(~type, ncol = 3, dir = "v") +
	geom_vline(xintercept=c(1), colour="white", size = 1.25) +
	geom_line(aes(x, y, color = fk), alpha = 0.9) +
	geom_vline(xintercept=c(1), colour="black", size = 0.75, linetype = "22") +
	scale_y_continuous(breaks = seq(0, 1, by=0.2), minor_breaks = seq(0, 1, by=0.1)) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
	coord_cartesian(xlim = c(-0.005, 2.7505), ylim = c(-0.005, 1.005), expand = F) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey90", "purple"))(3)[2],
																										 colorRampPalette(c("grey70", "navy"))(3)[2],
																										 colorRampPalette(c("grey20", "red"))(3)[2],
																										 colorRampPalette(c("grey90", "goldenrod3"))(3)[2])))(6)) +
	theme_few() +
	theme(aspect.ratio = 200/401,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(size = 1/3, colour = "grey80"),
				panel.grid.minor.y = element_line(size = 1/3, colour = "grey80"),
				panel.grid.major.x = element_line(size = 1/3, colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance to true breakpoint") +
	ylab("CDF") +
	guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
hist
ggsave(hist, filename = sprintf("_plot.hist.%s.pdf", prefix), width = 12, height = 7)





brk = seq(log10(1), log10(100e6), length.out = 81)

p = lapply(split(d, list(d$type, d$side)), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	#mat$fk  = x$fk[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


x = (p$side == "LHS")
p$Var1[x] = p$Var1[x] * -1


tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


tm = rbind(data.table(x = c(tik*-1 -1,tik+1), y = 0, xend = c(tik*-1 -1,tik+1), yend = 5),
					 data.table(y = tik, x = 80, yend = c(tik,tik), xend = 77),
					 data.table(y = tik, x = -80, yend = c(tik,tik), xend = -77))


scatter = ggplot(data = p) +
	facet_wrap(~type, ncol = 3, dir = "v") +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	geom_abline(slope = c(-1, 1), intercept = c(0,0), colour="black", alpha = 0.25) +
	#scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 2.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	#scale_fill_gradient2(low = "grey", mid = "darkblue", high = "orange", na.value = "grey85", midpoint = 1.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_fill_gradientn(colours = c("white", "snow", "yellow", "gold", "orange", "orangered1", "orangered3", "firebrick3", "firebrick4"), na.value = "grey85", limits=c(0,2.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 200/401,
				panel.background = element_rect(fill = "grey60"),
				legend.background = element_rect(fill = "grey85"),
				axis.text = element_text(size = 7),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				legend.title = element_blank()) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")
scatter
ggsave(scatter, filename = sprintf("_plot.scatter.%s.pdf", prefix), width = 12, height = 7)



### STATS 
se <- function(x) sqrt(var(x)/length(x))

rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}


x = by(d, list(d$type), function(x) cor(x$tru, x$det, method = "p")^2)
t(t(array(x, dim(x), dimnames(x))))

x = by(d, list(d$type), function(x) cor(x$tru, x$det, method = "s"))
t(t(array(x, dim(x), dimnames(x))))

x = by(d, list(d$type), function(x) rmsle(x$tru, x$det))
t(t(array(x, dim(x), dimnames(x))))



l = "R squared (Pearson)"
m = c(0,1)
x = by(d, list(d$fk, d$type), function(x) cor(x$tru, x$det, method = "p")^2)
t(array(round(x,4), dim(x), dimnames(x)))


l = "Spearman rank correlation coefficient"
m = c(0,1)
x = by(d, list(d$fk, d$type), function(x) cor(x$tru, x$det, method = "s"))
t(array(round(x,4), dim(x), dimnames(x)))


l = "RMSLE"
m = c(0,4)
x = by(d, list(d$fk, d$type), function(x) rmsle(x$tru, x$det))
t(array(round(x,4), dim(x), dimnames(x)))


ggplot(melt(array(x, dim(x), dimnames(x)))) +
	geom_line(aes(x = Var1, y = value, color = Var2)) +
	geom_point(aes(x = Var1, y = value, color = Var2), alpha = 0.5) +
	xlab("Focal allele count, fk") + ylab(l) +
	scale_x_continuous(breaks = c(2, 5, 10, 15, 20, 25)) +
	coord_cartesian(ylim = m) +
	theme_bw() +
	theme(legend.title = element_blank())
ggsave(filename = sprintf("_plot.stat_%s.%s.pdf", make.names(l), prefix), width = 8, height = 4)




x = 1

lhs = match$true.lhs.idx[x] + 1
rhs = match$true.rhs.idx[x] - 1

h0 = match$h0[x]
h1 = match$h1[x]

h0 = H[, h0]
h1 = H[, h1]

hh = rep("", length(POS) - 1)

for (i in 1:(length(POS) - 1)) {
	j = i+1
	hh[i] = sprintf("%d%d%d%d", h0[i], h0[j], h1[i], h1[j])
}

table(hh)

ins = (lhs):(rhs)
out = c(1:(lhs - 1), (rhs + 1):(length(POS)-1))

t = table(hh[ins]); round(t/sum(t), 4)

t = table(hh[out]); round(t/sum(t), 4)




### Length

DIST = c(0, cumsum(RATE[-(length(RATE))] * diff(POS) * 1e-6)) ### genetic distance
names(DIST) = as.character(POS)

h.gen = (DIST[ match(match$h.rhs, POS) ] - DIST[ match(match$h.lhs, POS) ]) + 1e-08
p.gen = (DIST[ match(match$p.rhs, POS) ] - DIST[ match(match$p.lhs, POS) ]) + 1e-08
g.gen = (DIST[ match(match$g.rhs, POS) ] - DIST[ match(match$g.lhs, POS) ]) + 1e-08

hh.gen = (DIST[ match(match$hh.rhs.pos, POS) ] - DIST[ match(match$hh.lhs.pos, POS) ]) + 1e-08
hp.gen = (DIST[ match(match$pp.rhs.pos, POS) ] - DIST[ match(match$pp.lhs.pos, POS) ]) + 1e-08
hg.gen = (DIST[ match(match$hmm.rhs, POS) ] - DIST[ match(match$hmm.lhs, POS) ]) + 1e-08

ch.gen = (DIST[ match(match$ch.rhs.pos, POS) ] - DIST[ match(match$ch.lhs.pos, POS) ]) + 1e-08
cc.gen = (DIST[ match(match$cc.rhs.pos, POS) ] - DIST[ match(match$cc.lhs.pos, POS) ]) + 1e-08

t.gen = (DIST[ match$true.rhs.idx ] - DIST[ match$true.lhs.idx ]) + 1e-08


q = rbind(data.table(idx = match$index, fk = match$fk, type = "True IBD",                    x = t.gen),
					data.table(idx = match$index, fk = match$fk, type = "(a) FGT, true haplotypes",    x = h.gen),
					data.table(idx = match$index, fk = match$fk, type = "(b) FGT, phased haplotypes",  x = p.gen),
					# data.table(idx = match$index, fk = match$fk, type = "(c) DGT, genotypes",          x = g.gen),
					data.table(idx = match$index, fk = match$fk, type = "(d) HMM, true haplotypes",    x = hh.gen),
					data.table(idx = match$index, fk = match$fk, type = "(e) HMM, phased haplotypes",  x = hp.gen),
					#data.table(idx = match$index, fk = match$fk, type = "(f) HMM, genotypes",          x = hg.gen),
					data.table(idx = match$index, fk = match$fk, type = "(g) HMM*, tru haplotypes",    x = ch.gen),
					data.table(idx = match$index, fk = match$fk, type = "(h) HMM*, phased haplotypes", x = cc.gen))

tl = unique(q$type)

q = split(q, list(q$type, q$fk))

q = lapply(q, function(x) {
	z = x
	w = fivenum(x$x)
	z$low = w[2]
	z$med = w[3]
	z$upp = w[4]
	z[1,]
})
q = rbindlist(q)


q$type = factor(q$type, levels = tl, ordered = T)

ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))

len = ggplot(q) + #[sample(1:nrow(p), 10000),]) +
	#facet_wrap(~mode, scales = "free", ncol = 1) +
	#geom_violin(aes(fk, x, fill = type), size = 0, alpha = 0.9) +
	#geom_pointrange(data=q, aes(x=fk, y=med, ymin=low, ymax=upp, group=type), color = "black", size = 2/3, shape='|', position=position_dodge(width = 0.9)) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=type), size = 1.5, position=position_dodge(width = 0.9)) +
	geom_point(aes(x=fk, y=med, group=type), size = 5, shape="-", position=position_dodge(width = 0.9)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), labels = format(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), scientific = F), minor_breaks = mb) +
	scale_x_continuous(breaks = 2:25) +
	#scale_color_manual(values = c("grey20", "royalblue1", "purple", "limegreen")) +
	#scale_fill_manual(values = c("turquoise3", "purple", "limegreen", "grey50")) +
	coord_cartesian(ylim = c(0.004, 12.5), xlim = c(1.5, 25.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "vertical",
				legend.text = element_text(size = 9),
				legend.key.height = unit(0.1, "cm"),
				legend.background = element_rect(fill = "white", colour = "grey50", size = 0.5),
				#panel.spacing.x = unit(-1, "points"), 
				strip.text = element_text(face = "bold", hjust = 0),
				#strip.text.y = element_blank(),
				#axis.ticks.y = element_blank(),
				panel.grid = element_line(colour = "grey80", size = 1/3),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") +
	ylab("Genetic length (cM)") +
	guides(color = guide_legend(override.aes = list(size = 5)))
len
ggsave(len, filename = sprintf("_plot.genlen.%s.pdf", prefix), width = 12, height = 9)






### phasing error

library(data.table)
library(ggplot2)
library(ggthemes)

load("data.OutOfAfricaHapMap20.H.RData")
load("data.OutOfAfricaHapMap20.P.RData")


x0 = 3
x1 = 4

xa = a = H[, x0]
xb = b = H[, x1]

het = which(a != b)

ff = het[as.integer(round(c(0.01, 0.011, 0.02, 0.03, 0.4, 0.5, 0.6, 0.97, 0.98, 0.99) * length(het)))]
ss = het[as.integer(round(c(0.1, 0.101, 0.45, 0.55, 0.9) * length(het)))]

a[ff] = xb[ff]
b[ff] = xa[ff]

for (s in ss) {
	rng = s:nrow(H)
	tmp = a[rng]
	a[rng] = b[rng]
	b[rng] = tmp
}

p0 = a
p1 = b

t0 = H[, x0]
t1 = H[, x1]


switch.flip.stat = function(p0, p1, t0, t1) {
	het = which(p0 != p1)
	n = length(p0)
	s = c()
	
	for (k in 1:(length(het) - 1)) {
		i = het[k]
		j = het[k+1]
		
		if (p0[i] != t0[i] && p0[j] != t0[j]) {
			s = c(s, i)
			
			rng = i:n
			tmp = p0[rng]
			p0[rng] = p1[rng]
			p1[rng] = tmp
		}
	}
	
	f = which(p0 != t0)
	
	list(s = s, f = f, ns = length(s), nf = length(f))
}


#switch.flip.stat(p0, p1, t0, t1)


prefix = "OutOfAfricaHapMap20"
#prefix = "OutOfAfricaHapMap20.GenErr_1000G"

load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.H.RData", prefix))
load(sprintf("data.%s.P.RData", prefix))


res = list()

for (i in 1:ncol(H)) {
	cat(i, "\n")
	if (i %% 2 == 0) next
	a = i
	b = i + 1
	
	p0 = P[, a]
	p1 = P[, b]
	t0 = H[, a]
	t1 = H[, b]
	
	res = c(res, list(switch.flip.stat(p0, p1, t0, t1)))
}



aac = table(AAC)

se <- function(x) sqrt(var(x)/length(x))


sw = rbindlist(lapply(res, function(x) {
	if (x$ns == 0) return(NULL)
	tmp = table(AAC[x$s])
	data.table(f = as.numeric(names(tmp)), p = as.vector(tmp) / as.vector(aac[names(tmp)]))
}))
sw = rbindlist(lapply(split(sw, sw$f), function(x) {
	data.table(f = x$f[1], mn = mean(x$p), se = se(x$p))
}))
plot(sw$f, sw$mn, log = 'y')


fl = rbindlist(lapply(res, function(x) {
	if (x$nf == 0) return(NULL)
	tmp = table(AAC[x$f])
	data.table(f = as.numeric(names(tmp)), p = as.vector(tmp) / as.vector(aac[names(tmp)]))
}))
fl = rbindlist(lapply(split(fl, fl$f), function(x) {
	data.table(f = x$f[1], mn = mean(x$p), se = se(x$p))
}))
plot(fl$f, fl$mn, log = 'y')



sw = lapply(res, function(x) {
	data.table(f = AAF[x$s], c = AAC[x$s])
})
sw = split(sw, sw$c)
sw = rbindlist(lapply(sw, function(x) {
	l = length(which(AAC == x$c[1]))
	n = nrow(x)
	data.table(f = x$f[1], c = x$c[1], p = n/l)
}))


ns = sapply(res, function(x) x$ns)
nf = sapply(res, function(x) x$nf)

s = unlist(lapply(res, function(x) x$s))
f = unlist(lapply(res, function(x) x$f))

ts = table(s)
ts = data.table(i = as.numeric(names(ts)), n = as.vector(ts))
ts$f = AAF[ts$i]
ts$c = AAC[ts$i]
ts = split(ts, ts$c)
ts = lapply(ts, function(x) {
	data.table(n = sum(x$n) , f = x$f[1])
})
ts = rbindlist(ts)

plot(ts$f, ts$n)


