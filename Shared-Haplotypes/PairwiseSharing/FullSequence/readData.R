
library("bigmemory")
options(bigmemory.typecast.warning=FALSE)


#args = commandArgs(TRUE)

markers.file = "f25.markers" #args[1]
sharing.file = "f25.sharing" #args[2]
genotypes.file = "f25.genotypes" #args[3]
panel.file = "f25.sample" #args[4]



# markers
M = read.table(markers.file, header = TRUE, stringsAsFactors = FALSE)


# sharing
tmp = read.table(sharing.file, header = TRUE, stringsAsFactors = FALSE)

tmp.as = rep('.', nrow(M))
tmp.fs = rep(0, nrow(M))
tmp.sh = rep(-1, nrow(M))
for (i in 1:nrow(tmp)) {
	x = tmp$ID[i] + 1
	tmp.as[x] = tmp$allele_shared[i]
	tmp.fs[x] = tmp$f[i]
	tmp.sh[x] = tmp$shared_haplotype[i]
}

M = cbind(M, allele_shared = tmp.as, f = tmp.fs, shared_haplotype = tmp.sh, stringsAsFactors = FALSE)

s = strsplit(tmp$sample_ids, "|", TRUE)
names(s) = tmp$ID + 1


# genotypes
I = strsplit(readLines(genotypes.file, n = 1), " ", TRUE)[[1]]

G = read.big.matrix(filename = genotypes.file, sep = " ", header = TRUE, col.names = I, ignore.row.names = TRUE, type = "char", 
										backingpath = getwd(),
										backingfile =  sprintf("R.%s.bin", genotypes.file), 
										descriptorfile = sprintf("R.%s.desc", genotypes.file))


S = big.matrix(nrow(G), ncol(G), type = "char", init = as.integer(0), dimnames = list(NULL, I),
							 backingpath = getwd(),
							 backingfile =  sprintf("R.%s.bin", sharing.file), 
							 descriptorfile = sprintf("R.%s.desc", sharing.file))

for (f in names(s)) {
	row = as.numeric(f)
	col = s[[f]]
	S[row, col] = as.integer(1)
}


# panel
P = read.table(panel.file, header = TRUE, stringsAsFactors = FALSE)


save(M, S, G, P, I, file = "R.f25.data.RData")


stop()


load("R.f25.data.RData")
G = attach.big.matrix(dget("R.f25.genotypes.desc"))
S = attach.big.matrix(dget("R.f25.sharing.desc"))

indv = colnames(G)




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


shared.sequence = function(g.target, g.other, site, break.max = 10) {
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
	list(shared=sh, breaks=br, abs.beg=abs.beg, abs.foc=abs.foc, abs.end=abs.end, rel.beg=rel.beg, rel.foc=rel.foc, rel.end=rel.end)
}



select.breakpoint.position = function(left, sh, m) {
	brks = sh$breaks
	brks[sh$rel.foc] = TRUE
	
	if (left) {
		brks = brks[(sh$rel.beg):(sh$rel.foc)]
		brks = rev(m$position[ which(brks) + sh$abs.beg - 1 ])
	} else {
		brks = brks[(sh$rel.foc):(sh$rel.end)]
		brks = m$position[ which(brks) + sh$abs.foc - 1 ]
	}

	med = median(brks)
	mad = mad(brks)
	dst = abs(brks - med) / mad
	brks = brks[ which(dst < mean(dst)) ]
	if (left) return(max(brks)) else return(min(brks))
	
	median.exclusion = function(brks) {
		l = length(brks)
		if (l == 2) {
			return(brks[1])
		}
		med = median(brks)
		mad = mad(brks)
		dist.beg = abs(brks[1] - med) / mad
		dist.end = abs(brks[l] - med) / mad
		
		if (dist.beg > dist.end) {
			return( median.exclusion(brks[-1]) )
		}
		if (dist.end > dist.beg) {
			return( median.exclusion(brks[-l]) )
		}
		return( brks[1] )
	}
	
	median.exclusion(brks)
}


make.search.matrix = function(beg, end, m, s, indv) {
	rng = beg:end
	rvs = rng[which(m$f[rng] != 0)]
	tag = as.character(rvs)
	mat = matrix(FALSE, nrow = length(rvs), ncol = length(indv), dimnames = list(tag, indv))
	tmp = lapply(s[tag], function(s, indv) { 
		x = rep(FALSE, length(indv)) 
		names(x) = indv
		x[s] = TRUE
		x
	}, indv)
 	for (rv in tag) {
 		mat[rv, ] = tmp[[rv]]
 	}
 	mat
}

search.matrix = function(tag, mat) {
	col = colnames(mat)
	fnd = list()
	for (row in rownames(mat)) {
		if (mat[row, tag]) {
			shr = col[as.logical(mat[row, ])]
			shr = shr[-(which(shr == tag))]
			if (length(shr) != 0) {
				fnd[[row]] = shr
			}
		}
	}
	fnd
}





site = sample(which(M$f !=0), 1)

s = S[[as.character(site)]]

target = sample(s, 1)
others = s[which(s != target)]


sharing = list()

for (other in others) {
	sharing[[other]] = shared.sequence(G[, target], G[, other], site)
}




library(ggplot2)

d = NULL
for (other in others) {
	hap = sharing[[other]]$shared
	nas = which(is.na(hap))
	if (length(nas) != 0) hap[nas] = -1
	
	hap[sharing[[other]]$breaks] = -2
	
	hap[sharing[[other]]$rel.foc] = 2
	
	d = rbind(d, data.frame(tag = other, pos = sharing[[other]]$abs.beg + (1:(sharing[[other]]$rel.end)) - 1, hap = as.character(hap), arg = hap / 6))
}

ggplot(data=d) + 
	geom_segment(aes(x=pos, xend=pos+1, y=tag, yend=as.numeric(tag) + arg, colour=hap)) 


d = NULL
for (other in others) {
	bl = select.breakpoint.position(TRUE,  sharing[[other]], M)
	br = select.breakpoint.position(FALSE, sharing[[other]], M)
	
	brks = M$position[ which(sharing[[other]]$breaks) + sharing[[other]]$abs.beg - 1 ]
	type = rep("0", length(brks))
	
	type[which(brks == bl)] = "1"
	type[which(brks == br)] = "1"
	
	d = rbind(d, data.frame(tag = other, 
													beg = M$position[ sharing[[other]]$abs.beg ], 
													foc = M$position[ sharing[[other]]$abs.foc ], 
													end = M$position[ sharing[[other]]$abs.end ], 
													brks = brks,
													type = type,
													stringsAsFactors = FALSE))
}

ggplot(data=d) + 
	geom_segment(aes(x=beg, xend=end, y=tag, yend=tag)) +
	geom_point(aes(x=foc, y=tag)) +
	geom_point(aes(x=brks, y=tag, alpha = d$type), shape="|", size = 5) +
	scale_alpha_manual(values=c("0"=0.2, "1"=1))





shr = list()

for (sub in others) {
	
	hsm = make.search.matrix(sharing[[sub]]$abs.beg, sharing[[sub]]$abs.end, M, S, indv)
	
	shr[[sub]] = search.matrix(sub, hsm)
	
}

if (length(shr) != 0) {
	csh = unlist(shr, recursive = TRUE, use.names = FALSE)
	sort(table(csh))
	ush = unique(csh)
}






target = sample(indv, 1)


shares = lapply(S, function(x, t) if (t %in% x) x , target)
shares = shares[ which(unlist(lapply(shares, length)) != 0) ]


for (rv in names(shares)) {
	cat(rv, "\n")
	
	others = shares[[rv]]
	shares[[rv]] = list()
	
	for (other in others) {
		if (other == target) {
			next
		}
		
		shares[[rv]][[other]] = shared.sequence(G[, target], G[, other], as.numeric(rv))
	}
	
	if (length(shares[[rv]]) == 0) {
		shares[[rv]] = NULL
	}
}




d = NULL
for (rv in names(shares)[1:100]) {
	for (other in names(shares[[rv]])) {
		sh = shares[[rv]][[other]]
		
		#bl = select.breakpoint.position(TRUE,  sh, M)
		#br = select.breakpoint.position(FALSE, sh, M)
		
		brks = M$position[ which(sh$breaks) + sh$abs.beg - 1 ]
		type = rep("0", length(brks))
		
		#type[which(brks == bl)] = "1"
		#type[which(brks == br)] = "1"
		
		d = rbind(d, data.frame(tag = other, 
														beg = M$position[ sh$abs.beg ], 
														foc = M$position[ sh$abs.foc ], 
														end = M$position[ sh$abs.end ], 
														brks = brks,
														type = type,
														stringsAsFactors = FALSE))
	}
}

library(ggplot2)

d = d[order(d$foc), ]

ggplot(data=d) + 
	geom_segment(aes(x=beg, xend=end, y=tag, yend=tag)) +
	geom_point(aes(x=foc, y=tag), shape="|", colour="black") +
	geom_point(aes(x=brks, y=tag), shape="|", colour = "white")





















