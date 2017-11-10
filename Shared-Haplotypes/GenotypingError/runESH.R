#
# Infer ESH
#


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
	if (g.target == 0) {
		if (g.other == 2) {
			return(-1)
		} else {
			return(0)
		}
	}
	if (g.target == 1) {
		if (g.other == 0) return(0)
		if (g.other == 1) return(NA)
		if (g.other == 2) return(1)
	}
	if (g.target == 2) {
		if (g.other == 0) {
			return(-1)
		} else {
			return(1)
		}
	}
	return(NULL)
}


chrom.shared.haplotype = function(g.target, g.other) {
	apply(matrix(c(g.target, g.other)), 1, shared.haplotype)
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
	
	br = which(br)
	
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


# completely random
gen.error.random = function(g, rate, M) {
	size = length(g)
	rate = round(size * rate)
	
	site = sample(1:size, size = rate)
	
	repl = sapply(g[site], function(x) {
		if (x == 0) return( x + 1 )
		if (x == 1) return( x + sample(c(-1, 1), 1) )
		if (x == 2) return( x - 1 )
	})
	
	list(site=site, repl=repl)
}

# frequency dependent
gen.error.frqdep = function(g, rate, M) {
	size = length(g)
	rate = round(size * rate)
	
	maf = apply(M[, c("af0", "af1")], 1, min)
	
	site = sample(1:size, size = rate, replace = FALSE, prob = 1 - maf)
	
	repl = sapply(g[site], function(x) {
		if (x == 0) return( x + 1 )
		if (x == 1) return( x + sample(c(-1, 1), 1) )
		if (x == 2) return( x - 1 )
	})
	
	list(site=site, repl=repl)
}

# frequency dep., but with homoz. and heteroz. separated
gen.error.homhet = function(g, rate, M) {
	size = length(g)
	rate = round(size * rate)
	
	maf = apply(M[, c("af0", "af1")], 1, min)
	
	m = g
	x = which(M$af0 < M$af1)
	if (length(x) != 0) {
		m[x] = sapply(m[x], function(m) {
			if (m == 0) return(2)
			if (m == 2) return(0)
			return(1)
		})
	}
	
	b0 = (m == 0)
	b1 = (m == 1)
	b2 = (m == 2)
	
	r0 = length(which(b0)) / size
	r1 = length(which(b1)) / size
	r2 = length(which(b2)) / size
	
	prob = rep(0, size)
	prob[b0] = qexp(1-maf[b0], r0) # major hom
	prob[b1] = qexp(1-maf[b1], r1) # het
	prob[b2] = qexp(1-maf[b2], r2) # minor hom
	
	zygotype = c("0"="major hom", "1"="heterozyg", "2"="minor hom")
	d = data.frame(weight = prob / sum(prob), maf = maf, type = zygotype[as.character(m)])
	ggplot(data = d) + geom_line(aes(x=maf, y=weight, colour=type))
	
	site = unique(sample(1:size, size = rate, replace = TRUE, prob = prob))
	
	repl = sapply(g[site], function(x) {
		if (x == 0) return( x + 1 )
		if (x == 1) return( x + sample(c(-1, 1), 1) )
		if (x == 2) return( x - 1 )
	})
	
	list(site=site, repl=repl)
}


prepare.error = function(gen.file, error.rates, error.funcs, M) {
	G = load.bigmatrix(gen.file)
	
	tags = colnames(G)
	
	err = list()
	
	for (tag in tags) {
		cat(tag, "\n")
		err[[tag]] = list()
		
		for (func in names(error.funcs)) {
			
			err[[tag]][[func]] = list()
			error.func = error.funcs[[func]]
			
			for (error.rate in error.rates) {
				
				rate = as.character(error.rate * 100)
				err[[tag]][[func]][[rate]] = error.func(g, error.rate, M)
				
			}
		}
	}
	
	err
}


apply.error = function(g, indv, func, rate, error) {
	err = error[[indv]][[func]][[rate]]
	g[err$site] = err$repl
	g
}





target.csh = function(target, shared, error, M, gen.file) {
	
	
	list()
}






esh.segments = function(target, shared, M, gen.file, error.rate, error.func, mc.cores = NULL) {
	scan = function(rvs, target, shared, M, gen.file, error.rate, error.func) {
		G = load.bigmatrix(gen.file)
		
		rv = as.numeric(rvs)
		others = shared[[rvs]]
		
		g.target = G[, target]
		
		if (error.rate != 0) {
			g.target = error.func(g.target, error.rate, M)
		}
		
		esh = list()
		
		for (other in others) {
			
			g.other = G[, other]
			
			if (error.rate != 0) {
				g.other = error.func(g.other,  error.rate, M)
			} 
			
			# check if rare variant is still detectable
			if (g.target[rv] == 1 && g.other[rv] == 1) { 
				esh[[other]] = shared.sequence(g.target, g.other, rv, M, break.max = nrow(M))
			} else {
				esh[[other]] = NA
			}
			
		}
		
		esh
	}
	
	esh = list()
	
	if (is.null(mc.cores) || is.na(mc.cores)) {
		
		i = 1
		l = length(shared)
		
		for (rvs in names(shared)) {
			if (i %% 100 == 0) {
				cat(sprintf(" %.1f%%\n", (i/l)*100))
			}
			i = i + 1
			
			esh[[rvs]] = scan(rvs, target, shared, M, gen.file, error.rate, error.func)
		}
		
	} else {
		rvs = names(shared)
		names(rvs) = names(shared)
		
		esh = mclapply(rvs, scan, target, shared, M, gen.file, error.rate, error.func, mc.cores = mc.cores)
	}
	
	esh
}



###


library(bigmemory)
library(parallel)



# data.file = "data.sim.2000.RData"
# save.dir = "result_ESH"
# cores = NULL


args = commandArgs(TRUE)
data.file = args[1]
save.dir = args[2]
cores = args[3]


dir.create(file.path(getwd(), save.dir), showWarnings = FALSE)


load(data.file)

#cat("Determine sharing ... ")
#sharing = sharing.list(I, shh.file)
#cat("OK\n")

load("data.sharing.RData")


# error settings
error.rates = c(0, exp(seq(log(0.01), log(10), length.out = 20))) / 100
error.funcs = list("random"=gen.error.random, "frqdep"=gen.error.frqdep, "homhet"=gen.error.homhet)
# !


i = 0
l = length(I)
for (target in sample(I[1:100])) {
	i = i + 1
	
	for (error.rate in sample(error.rates)) {
		
		rate.tag = as.character(error.rate * 100)
		
		for (func.tag in sample(names(error.funcs))) {
			
			cat(i, "of", l, "-->", J[target], func.tag, rate.tag, "\n")
			
			save.file = file.path(getwd(), sprintf("%s/result.ESH.%s__%s__%s.RData", save.dir, func.tag, rate.tag, J[target]))
			
			if (file.exists(save.file)) {
				cat("EXISTS\n")
				next
			}
			
			ESH = esh.segments(target, sharing[[target]], M, gen.file, error.rate, error.funcs[[func.tag]], mc.cores = cores)
			
			if (file.exists(save.file)) {
				cat("EXISTS\n")
				next
			}
			
			save(ESH, file = save.file)
			
			cat("\n")
			
		}
	}
	
}







