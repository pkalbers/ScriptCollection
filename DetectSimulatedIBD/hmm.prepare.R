#
# get empirical HMM probs
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
library(ggplot2)
library(ggthemes)


args = commandArgs(T)

prefix = args[1] # "history.generror_1000g"

tru.file = args[2] # "../truth/ibd_mrca.truth.RData"

load(tru.file)


del = union(which(truth$wall.lhs), which(truth$wall.rhs))
if (length(del) > 0) {
	truth = truth[-del, ]
}


load(sprintf("%s.RData", prefix))


del = which(FKI$frq > 0.5)
if (length(del) > 0) {
	del = FKI$index[del]
	del = which(truth$index %in% del)
	if (length(del) > 0) {
		truth = truth[-del, ]
	}
}

G = as.matrix(load.bigmatrix(sprintf("%s.G", prefix)))

AC = rowSums(G)
AF = AC / (ncol(G) * 2)



G0 = rep(0, nrow(truth))
i = which(truth$h0 %% 2 != 0)
j = which(truth$h0 %% 2 == 0)
G0[i] = (truth$h0[i] + 1) / 2
G0[j] = truth$h0[j] / 2

G1 = rep(0, nrow(truth))
i = which(truth$h1 %% 2 != 0)
j = which(truth$h1 %% 2 == 0)
G1[i] = (truth$h1[i] + 1) / 2
G1[j] = truth$h1[j] / 2



### inital state probs

FS = truth$index # focal sites
FK = AC[FS] # focal freq/k

c0 = as.matrix(data.frame(FS, G0))
c1 = as.matrix(data.frame(FS, G1))

g0 = G[c0]
g1 = G[c1]

d = data.table(fk = FK, ibd = (g0 == 1 & g1 == 1))
d = split(d, d$fk)
d = lapply(d, function(x) {
	z = x[1, ]
	z$ibd = length(which(x$ibd)) / nrow(x)
	z
})
d = rbindlist(d)

del = which(d$fk < 2 | d$fk > max(FKI$n.sharer))
if (length(del) > 0) {
	d = d[-del, ]
}


gg = ggplot(d) + 
	geom_line(aes(fk, ibd)) + 
	geom_point(aes(fk, ibd)) + 
	scale_x_continuous(breaks = 2:max(FKI$n.sharer)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
	coord_cartesian(xlim = c(2, max(FKI$n.sharer)), ylim = c(0, 1)) +
	theme_bw() +
	xlab("Rare allele count (fk)") +
	ylab("Proportion of correctly observed focal genotype pairs")

ggsave(gg, filename = sprintf("_plot.hmm.inital.%s.pdf", prefix), width = 10, height = 10)

hmm.init = d



### emission probs

max = length(POS)

ins = matrix(0, nrow = max, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
out = ins

trans = rep(NA, 23)
trans[c(1:3, 11:13, 21:23)] = c(1, 2, 3,
																2, 4, 5,
																3, 5, 6)

n = 0
lim = min(5000, nrow(truth))

n.indv = ncol(G)

cat("Empirical emissions, sampled ...\n")

for (i in sample(nrow(truth))[1:lim]) {
	n = n + 1
	if (n %% 1000 == 0) cat(sprintf(" %d of %d\n", n, lim))
	
	lhs = match(truth$lhs[i], POS)
	rhs = match(truth$rhs[i], POS)
	foc = truth$index[i]
	
	rng = c(lhs : (foc - 1), (foc + 1) : rhs)
	sgm = lhs : rhs
	
	s0 = G0[i]
	s1 = G1[i]
	
	share.pair = trans[ ((G[rng, s0] * 10) + (G[rng, s1] + 1)) ]
	
	rares = sgm[ which(AC[sgm] >= 1 & AC[sgm] <= 25) ]
	o0 = sample(n.indv, size = 1)
	o1 = sample(n.indv, size = 1)
	while(T) {
		if (o0 != o1) {
			tmp = which(G[rares, o0] == 1 & G[rares, o1] == 1)
			if (length(tmp) == 0) {
				break
			}
		}
		o0 = sample(n.indv, size = 1)
		o1 = sample(n.indv, size = 1)
	}
	
	other.pair = trans[ ((G[rng, o0] * 10) + (G[rng, o1] + 1)) ]
	
	coord = matrix(c(rng, share.pair), ncol = 2, byrow = F)
	ins[ coord ] = ins[ coord ] + 1
	
	coord = matrix(c(rng, other.pair), ncol = 2, byrow = F)
	out[ coord ] = out[ coord ] + 1
	
	
# 	g0 = G0[i]
# 	g1 = G1[i]
# 	
# 	g = trans[ ((G[, g0] * 10) + (G[, g1] + 1)) ]
# 	
# 	lhs = match(truth$lhs[i], POS)
# 	rhs = match(truth$rhs[i], POS)
# 	foc = truth$index[i]
# 	
# 	r.ins = c(lhs : (foc - 1), (foc + 1) : rhs)
# 	r.out = c(  1 : (lhs - 1), (rhs + 1) : max)
# 	
# 	coord = matrix(c(r.ins, g[r.ins]), ncol = 2, byrow = F)
# 	ins[ coord ] = ins[ coord ] + 1
# 	
# 	coord = matrix(c(r.out, g[r.out]), ncol = 2, byrow = F)
# 	out[ coord ] = out[ coord ] + 1
}

#save(ins, out, file = sprintf("_hmm.emission.counts.%s.RData", prefix))
#load(sprintf("_hmm.emission.counts.%s.RData", prefix))

d.ins = matrix(0, nrow = length(unique(AC)), ncol = 6, dimnames = list(as.character(sort(unique(AC))), c("00", "01", "02", "11", "12", "22")))
d.out = d.ins

for (ac in sort(unique(AC))) {
	over = which(AC == ac)
	
	tmp = if (length(over) == 1) ins[over, ] else colSums(ins[over, ])
	tmp = tmp / sum(tmp)
	d.ins[as.character(ac), ] = tmp
	
	tmp = if (length(over) == 1) out[over, ] else colSums(out[over, ])
	tmp = tmp / sum(tmp)
	d.out[as.character(ac), ] = tmp
}


# d = rbind(cbind(type =  "inside", melt(d.ins, as.is = T)), 
# 					cbind(type = "outside", melt(d.out, as.is = T)))
# 
# ggplot(data = d) + 
# 	facet_grid(~type) + 
# 	geom_raster(aes(x = factor(Var2), y = as.numeric(Var1) / 5000, fill = value)) +
# 	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 0.5) +
# 	coord_cartesian(expand = F) +
# 	theme_base()
# 
# ggplot(data = d) + 
# 	facet_grid(~Var2) + 
# 	geom_raster(aes(x = type, y = as.numeric(Var1) / 5000, fill = value)) +
# 	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 0.6) +
# 	coord_cartesian(expand = F) +
# 	theme_base()


d.ins -> xxx
d.ins <- xxx

d.ins[, '00'] = d.ins[, '00'] * (4/4)
d.ins[, '01'] = d.ins[, '01'] * (3/4)
d.ins[, '02'] = d.ins[, '02'] * (0.01/4)
d.ins[, '11'] = d.ins[, '11'] * (2/4)
d.ins[, '12'] = d.ins[, '12'] * (3/4)
d.ins[, '22'] = d.ins[, '22'] * (4/4)

range(rowSums(d.ins))

d = rbind(cbind(melt(d.out, as.is = T), type = "outside"),
					cbind(melt(d.ins, as.is = T), type =  "inside"))
names(d) = c("freq", "pair", "prop", "type")
d$freq = as.numeric(d$freq) / 5000


q = sort(unique(d$freq)) # (0:5000)/5000
p = 1-q

e.out = rbind(data.table(freq = q, pair = "00", prop = p^4),
							data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
							data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
							data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
							data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
							data.table(freq = q, pair = "22", prop = q^4))
e.out$type = "outside"

e.ins = rbind(data.table(freq = q, pair = "00", prop = p^3),
							data.table(freq = q, pair = "01", prop = 2*(p^2)*q),
							data.table(freq = q, pair = "02", prop = 0),
							data.table(freq = q, pair = "11", prop = ((p^2)*q) + (p*(q^2))),
							data.table(freq = q, pair = "12", prop = 2*p*(q^2)),
							data.table(freq = q, pair = "22", prop = q^3))
e.ins$type = "inside"

e = rbind(e.out, e.ins)

x = d
x$prop = d$prop - e$prop

gg = ggplot(data = d) + 
	facet_grid(~type) + 
	#geom_smooth(aes(y = prop, x = freq, colour = pair)) +
	geom_line(aes(y = prop, x = freq, colour = pair)) +
	geom_line(data = e, aes(y = prop, x = freq, group  = pair), alpha = 0.75, colour = "black", linetype="dashed") +
	theme_gdocs() +
	theme(legend.title = element_blank()) +
	xlab("Allele frequency") +
	ylab("Observed proportions of genotype pairs")

ggsave(gg, filename = sprintf("_plot.hmm.emission.%s.pdf", prefix), width = 10, height = 10)


hmm.emiss.ibd = d.ins
hmm.emiss.non = d.out


save(hmm.init, hmm.emiss.ibd, hmm.emiss.non, file = sprintf("hmm.probs.%s.RData", prefix))






### emission probs ?????

hdst.ibd = NULL
hdst.non = NULL

n = 0
lim = min(5000, nrow(truth))

n.indv = ncol(G)

cat("Empirical emissions, sampled ...\n")

for (i in sample(nrow(truth))[1:lim]) {
	n = n + 1
	if (n %% 1000 == 0) cat(sprintf(" %d of %d\n", n, lim))
	
	lhs = match(truth$lhs[i], POS)
	rhs = match(truth$rhs[i], POS)
	foc = truth$index[i]
	
	rng = c(lhs : (foc - 1), (foc + 1) : rhs)
	sgm = lhs : rhs
	
	s0 = G0[i]
	s1 = G1[i]
	
	share.pair = sum(abs(G[rng, s0] - G[rng, s1]))
	
	rares = sgm[ which(AC[sgm] >= 1 & AC[sgm] <= 25) ]
	o0 = sample(n.indv, size = 1)
	o1 = sample(n.indv, size = 1)
	while(T) {
		if (o0 != o1) {
			tmp = which(G[rares, o0] == 1 & G[rares, o1] == 1)
			if (length(tmp) == 0) {
				break
			}
		}
		o0 = sample(n.indv, size = 1)
		o1 = sample(n.indv, size = 1)
	}
	
	other.pair = sum(abs(G[rng, o0] - G[rng, o1]))
	
	hdst.ibd = rbind(hdst.ibd, data.table(fk = truth$fk[i], size = POS[rhs] - POS[lhs], hdst = share.pair, hnum = length(rng)))
	hdst.non = rbind(hdst.non, data.table(fk = truth$fk[i], size = POS[rhs] - POS[lhs], hdst = other.pair, hnum = length(rng)))
}

hdst.ibd$state = "IBD"
hdst.non$state = "NON"

d = rbind(hdst.ibd, hdst.non)
d$fk = factor(d$fk)

ggplot(d) + geom_point(aes(x = fk, y = hdst/hnum, colour = state), shape=1)


ggplot(d) + geom_point(aes(x = fk, y = hdst/size, colour = state), position = position_dodge(width = 0.5))
















