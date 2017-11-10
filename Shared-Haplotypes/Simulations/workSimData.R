

library(bigmemory)
library(parallel)
library(ggplot2)
library(grid)



sharing.list = function(I, shh.file) {
	S = load.bigmatrix(shh.file)
	sharing = list()
	for (target in I) {
		sharing[[target]] = list()
		for (rv in which(S[, target] == 1)) {
			others = names(which(S[rv, ] == 1))
			others = others[-(which(others == target))]
			sharing[[target]][[ as.character(rv) ]] = others
		}
	}
	sharing
}



shared.haplotype = function(g.target, g.other) {
	h = NA
	if (g.target == 0) {
		#if (g.other == 0) h = 0
		#if (g.other == 1) h = 0
		#if (g.other == 2) h = 0
		h = 0
	}
	if (g.target == 1) {
		if (g.other == 0) h = 0
		if (g.other == 1) h = NA
		if (g.other == 2) h = 1
	}
	if (g.target == 2) {
		#if (g.other == 0) h = 1
		#if (g.other == 1) h = 1
		#if (g.other == 2) h = 1
		h = 1
	}
	h
}



shared.breakpoint = function(g.target, g.other) {
	((g.target == 0 && g.other == 2) || (g.target == 2 && g.other == 0))
}



shared.sequence.side = function(left, g.target, g.other, site, break.max = 10) {
	sites = if(left) { site : 1 } else { site : length(g.target) }
	l = length(sites)
	n = 0
	i = 0
	k = 0
	sh = rep(NA, l)
	br = rep(NA, l)
	for (i in 1:l) {
		k = sites[i]
		
		if (i == 1) {
			if (g.other[k] == 1 || g.target[k] == 1) {
				sh[i] = 1
				br[i] = FALSE
				next
			}
		}
		
		sh[i] = shared.haplotype(g.target[k], g.other[k])
		br[i] = shared.breakpoint(g.target[k], g.other[k])
		
		if (br[i]) {
			n = n + 1
			if (n == break.max + 1) {
				break
			}
		}
	}
	
	sh = sh[1:i]
	br = br[1:i]
	
	list(shared = sh, breaks = br, site = k)
}



shared.sequence = function(g.target, g.other, site, M, break.max = 10) {
	l = shared.sequence.side(left = TRUE,  g.target, g.other, site, break.max)
	r = shared.sequence.side(left = FALSE, g.target, g.other, site, break.max)
	
	sh = c(rev(l$shared[-1]), r$shared)
	br = c(rev(l$breaks[-1]), r$breaks)
	
	abs.beg = l$site
	abs.foc = site
	abs.end = r$site
	
	rel.beg = 1
	rel.foc = length(l$shared)
	rel.end = length(sh)
	
	pos.beg = M$position[ l$site ]
	pos.foc = M$position[ site ]
	pos.end = M$position[ r$site ]
	
	list(shared=sh, breaks=br, 
			 abs.beg=abs.beg, abs.foc=abs.foc, abs.end=abs.end, 
			 rel.beg=rel.beg, rel.foc=rel.foc, rel.end=rel.end,
			 pos.beg=pos.beg, pos.foc=pos.foc, pos.end=pos.end)
}



esh.segments = function(target, shared, M, gen.file, mc.cores = NULL) {
	scan = function(rvs, target, shared, M, gen.file) {
		G = load.bigmatrix(gen.file)
		
		rv = as.numeric(rvs)
		others = shared[[rvs]]
		
		esh = list()
		for (other in others) {
			esh[[other]] = shared.sequence(G[, target], G[, other], rv, M, break.max = 10)
		}
		
		return(esh)
	}
	
	esh = list()
	
	if (is.null(mc.cores)) {
		for (rvs in names(shared)) {
			esh[[rvs]] = scan(rvs, target, shared, M, gen.file)
		}
	} else {
		rvs = names(shared)
		names(rvs) = names(shared)
		esh = mclapply(rvs, scan, target, shared, M, gen.file, mc.cores = mc.cores)
	}
	
	esh
}



ibd.segments = function(target, shared, M, hap.file, ibd.file, mc.cores = NULL) {
	scan = function(ibd, M, hap.file, ibd.file) {
		H = load.bigmatrix(hap.file)
		IBD = load.bigmatrix(ibd.file)
		
		rv = as.numeric(ibd$rvs)
		
		pair = strsplit(ibd$other, " ", TRUE)[[1]]
		
		is.mat = (H[rv, pair[1]] == 1)
		is.pat = (H[rv, pair[2]] == 1)
		
		if (is.mat && is.pat) {
			return(ibd)
		}
		
		ibd$other.hap = ifelse(is.mat, pair[1], pair[2])
		ibd$other.lab = ifelse(is.mat, "M", "P")
		
		lookup = as.list(sort(c(as.numeric(ibd$target.hap), as.numeric(ibd$other.hap))))
		match = mwhich(IBD, c("A","B"), lookup, list("eq", "eq"), "AND")
		
		if (length(match) == 0) {
			return(ibd)
		}
		
		sub = IBD[match, ]
		
		ibd$foc = M$position[rv]
		match = which(sub[, "beg"] <= ibd$foc & sub[, "end"] >= ibd$foc)
		
		if (length(match) == 0) {
			return(ibd)
		}
		
		if (length(match) > 1) {
			sub = sub[match, ]
			dst = sub[, "end"] - sub[, "beg"]
			match = which.min(dst)
		} 
		
		ibd$beg = sub[match, "beg"]
		ibd$end = sub[match, "end"]
		
		ibd
	}
	
	H = load.bigmatrix(hap.file)
	
	ibd = NULL
	
	for (rvs in names(shared)) {
		rv = as.numeric(rvs)
		
		pair = strsplit(target, " ", TRUE)[[1]]
		
		is.mat = (H[rv, pair[1]] == 1)
		is.pat = (H[rv, pair[2]] == 1)
		
		if (is.mat && is.pat) {
			next
		}
		
		target.hap = ifelse(is.mat, pair[1], pair[2])
		target.lab = ifelse(is.mat, "M", "P")
		
		ibd = rbind(ibd, data.frame(rvs = rvs,
																target = target, other = shared[[rvs]], 
																target.hap = target.hap, other.hap  = NA, 
																target.lab = target.lab, other.lab  = NA, 
																beg = NA, foc = NA, end = NA))
	}
	
	ibd = apply(ibd, 1, as.list)

	if (is.null(mc.cores)) {
		for (i in 1:length(ibd)) {
			ibd[[i]] = scan(ibd, M, hap.file, ibd.file)
		}
	} else {
		ibd = mclapply(ibd, scan, M, hap.file, ibd.file, mc.cores = mc.cores)
	}
	
	out = NULL
	
	for (i in 1:length(ibd)) { 
		out = rbind(out, as.data.frame(ibd[[i]], stringsAsFactors = FALSE))
	}
	
	out
}


old.ibd.segments = function(seg, target, M, L, hap.file, ibd.file) {
	H = load.bigmatrix(hap.file)
	IBD = load.bigmatrix(ibd.file)
	
	rv = as.numeric(seg$rvs[1])
	ht = strsplit(target, " ", TRUE)[[1]]
	
	is.mat = as.logical(H[rv, ht[1]])
	is.pat = as.logical(H[rv, ht[2]])
	
	if (is.mat && is.pat) stop("Double in target")
	if (is.mat) ht = ht[1]
	if (is.pat) ht = ht[2]
	
	ibd = data.frame(target = as.character(target), ab = ifelse(is.mat, "M", "P"), ht = ht, other = unique(seg$tag), ho = NA, beg = NA, end = NA, stringsAsFactors = FALSE)
	
	for (i in 1:nrow(ibd)) {
		ho = strsplit(ibd$other[i], " ", TRUE)[[1]]
		
		is.mat = as.logical(H[rv, ho[1]])
		is.pat = as.logical(H[rv, ho[2]])
		
		if (is.mat && is.pat) stop("Double in other")
		if (is.mat) ibd$ho[i] = ho[1]
		if (is.pat) ibd$ho[i] = ho[2]
	}
	
	if (any(is.na(ibd$ho))) {
		print(ibd)
		ibd = ibd[-(which(is.na(ibd$ho))), ]
	}
	
	for (i in 1:nrow(ibd)) {
		ha = min(as.numeric(ibd$ht[i]), as.numeric(ibd$ho[i]))
		hb = max(as.numeric(ibd$ht[i]), as.numeric(ibd$ho[i]))
		ab = mwhich(IBD, c("A","B"), list(ha, hb), list("eq", "eq"), "AND")
		
		if (length(ab) == 0) next
		
		rng = IBD[ab, ]
		
		rvp = M$position[rv] * L
		k = which(rng[, "beg"] <= rvp & rng[, "end"] >= rvp)
		
		if (length(k) != 1) { 
			print(k)
			next 
		}
		
		ibd$beg[i] = rng[k, "beg"] / L
		ibd$end[i] = rng[k, "end"] / L
	}
	
	if (any(is.na(ibd$beg))) {
		print(ibd)
		ibd = ibd[-(which(is.na(ibd$beg))), ]
	}
	
	list(seg=seg, ibd=ibd)
}



###


args = commandArgs(TRUE)
file = args[1]


load(sprintf("data.%s.RData", file))


sharing = sharing.list(I, shh.file)

result = list()

for (target in sample(I, 5)) {
	cat(target, "\n")
	
	shared = sharing[[target]]
	
	esh = esh.segments(target, shared, M, gen.file, mc.cores = 16)
	ibd = ibd.segments(target, shared, M, hap.file, ibd.file, mc.cores = 16)
	
	result[[target]] = list(esh=esh, ibd=ibd)
}

save(result, file = sprintf("result.%s.RData", file))




stop()




d = NULL
for (tag in names(ESH)) {
	for (other in names(ESH[[tag]])) {
		esh = ESH[[tag]][[other]]
		brk = NA
		if (any(esh$breaks)) {
			brk = M$position[ which(esh$breaks) + esh$abs.beg - 1 ]
		}
		d = rbind(d, data.frame(rvs = tag,
														tag = other, 
														beg = M$position[ esh$abs.beg ], 
														foc = M$position[ esh$abs.foc ], 
														end = M$position[ esh$abs.end ], 
														brk = brk,
														stringsAsFactors = FALSE))
	}
}



splt = split(d, d$rvs)

data = mclapply(splt, ibd.segment, target, M, L, hap.file, ibd.file, mc.cores = 6)


d = NULL
for (i in 1:length(data)) {
	seg = data[[i]]$seg
	ibd = data[[i]]$ibd
	
	for (k in 1:nrow(ibd)) {
		m = which(seg$tag == ibd$other[k])
		if (length(m) == 0) next
		d = rbind(d, data.frame(rvs = as.numeric(seg$rvs[m]), 
														focus = seg$foc[m],
														target = ibd$target[k],
														other = ibd$other[k],
														seg.beg = seg$beg[m],
														seg.end = seg$end[m],
														breaks = seg$brk[m],
														ibd.beg = ibd$beg[k],
														ibd.end = ibd$end[k],
														ibd.hap = ibd$ab[k],
														stringsAsFactors = FALSE))
	}
}
d = d[order(d$rvs), ]

d = rbind(data.frame(rvs = 0, focus = unique(d$focus), target=NA, other=as.character(target), seg.beg=0, seg.end=0, breaks=0, ibd.beg=0, ibd.end=0, ibd.hap="A", stringsAsFactors = FALSE), d)


gg = ggplot(data=d) + 
	facet_grid(rvs~., scales = "free_y", space = "free_y") +
	geom_segment(data = d, aes(x=ibd.beg, xend=ibd.end, y=other, yend=other, colour=ibd.hap), size=3) +
	geom_segment(aes(x=seg.beg, xend=seg.end, y=other, yend=other)) +
	geom_point(aes(x = -0.015, y=other, fill=other), colour="darkgrey", shape=22, size=4) +
	geom_point(aes(x=focus, y=other), size=2) +
	geom_point(aes(x=breaks, y=other), shape=124, size=3) +
	geom_vline(xintercept=c(0,1)) +
	theme_classic() +
	theme(axis.ticks.length=unit(0.5, "lines"),
				axis.title.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.line.y=element_blank(),
				strip.text.y=element_blank(),
				strip.background=element_blank(),
				panel.margin.x=unit(0.5, "lines"),
				legend.position="none") +
	scale_x_continuous(limits=c(-0.025, 1.025), expand = c(0,0), breaks=seq(0, 1, length.out = 11), labels=seq(0, L, length.out = 11) / 1e06) +
	scale_y_discrete(labels=J) +
	scale_fill_manual(values=K) +
	xlab("Position (Mb)")


ggsave(gg, filename = sprintf("plot.%s.png", names(target)), width = 10, height = length(unique(paste(d$rvs, d$focus, d$other))) / 7, limitsize=FALSE)




stop()



d = NULL
for (i in 1:length(data)) {
	ibd = data[[i]]$ibd
	for (k in 1:nrow(ibd)) {
		d = rbind(d, data.frame(target = ibd$target[k],
														other = ibd$other[k],
														beg = round(ibd$beg[k] * L),
														end = round(ibd$end[k] * L),
														hap = ibd$ab[k],
														stringsAsFactors = FALSE))
	}
}
d = unique(d)

profile.a = rep(0, L)
profile.b = rep(0, L)
for (i in 1:nrow(d)) {
	rng = (d$beg[i]):(d$end[i])
	if (d$hap[i] == "A") profile.a[rng] = profile.a[rng] + 1
	if (d$hap[i] == "B") profile.b[rng] = profile.b[rng] + 1
}

profile.a = data.frame(x = 1:L, y = profile.a, ab = rep("A", L))
profile.b = data.frame(x = 1:L, y = profile.b, ab = rep("B", L))

rng0 = 1:(nrow(profile.a)-1)
rng1 = 2:nrow(profile.a)
k = which(profile.a$y[rng0] != profile.a$y[rng1])
k = sort(c(1, L, k, k-1))
profile.a = profile.a[k, ]

rng0 = 1:(nrow(profile.b)-1)
rng1 = 2:nrow(profile.b)
k = which(profile.b$y[rng0] != profile.b$y[rng1])
k = sort(c(1, L, k, k-1))
profile.b = profile.b[k, ]

d = rbind(profile.a, profile.b)

ggplot(data=d) + geom_line(aes(x = x, y = y, colour = ab))
















