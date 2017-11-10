#
# analyse all results
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

#args = commandArgs(T)

prefix = "history"  #args[1]             # file prefix

load(sprintf("%s.RData", prefix))



end.lhs = round(POS[1], 1)
end.rhs = round(POS[length(POS)], 1)




res.files = dir(pattern = "^result\\..+\\.RData$", recursive = T)

data = NULL

for (res.file in res.files) {
	load(res.file)
	
	pack = sub("^.+\\.([0-9]+)\\..+$", "\\1", res.file)
	
	origin = sub("^([^\\/]*)\\/?pairs\\/.+", "\\1", res.file)
	if (origin == "") origin = "simulated"
	pair$origin = origin
	
	cat(sprintf("%s :: %s", origin, pack))
	
	
	fail.g = which(pair$g.lhs == 0 & pair$g.rhs == 0)
	fail.h = which(pair$h.lhs == 0 & pair$h.rhs == 0)
	fail.p = which(pair$p.lhs == 0 & pair$p.rhs == 0)
	fail = unique(c(fail.g, fail.h, fail.p))
	if (length(fail) > 0) {
		cat(sprintf(" [Failed:  %.3f%%]", length(fail) / nrow(pair) * 100))
		pair = pair[-fail, ]
		if (nrow(pair) == 0) {
			next
		}
	}
	
	
	pair$g.end.lhs = F; end = which(round(pair$g.lhs, 1) == end.lhs); if (length(end) > 0) pair$g.end.lhs[end] = T
	pair$g.end.rhs = F; end = which(round(pair$g.rhs, 1) == end.rhs); if (length(end) > 0) pair$g.end.rhs[end] = T
	
	pair$h.end.lhs = F; end = which(round(pair$h.lhs, 1) == end.lhs); if (length(end) > 0) pair$h.end.lhs[end] = T
	pair$h.end.rhs = F; end = which(round(pair$h.rhs, 1) == end.rhs); if (length(end) > 0) pair$h.end.rhs[end] = T
	
	pair$p.end.lhs = F; end = which(round(pair$p.lhs, 1) == end.lhs); if (length(end) > 0) pair$p.end.lhs[end] = T
	pair$p.end.rhs = F; end = which(round(pair$p.rhs, 1) == end.rhs); if (length(end) > 0) pair$p.end.rhs[end] = T
	
	
	pair$rel.g.lhs = pair$g.lhs - pair$position
	pair$rel.g.rhs = pair$g.rhs - pair$position
	
	pair$rel.h.lhs = pair$h.lhs - pair$position
	pair$rel.h.rhs = pair$h.rhs - pair$position
	
	pair$rel.p.lhs = pair$p.lhs - pair$position
	pair$rel.p.rhs = pair$p.rhs - pair$position
	
	
	cat("\n")
	
	data = rbind(data, pair)
}


###########


library(ggplot2)


get.gen.states = function(brk, g.pair0, g.pair1, S, E) {
	coord = as.matrix(data.table(x = c(brk, brk), y = c(g.pair0, g.pair1)))
	colnames(coord) = NULL
	
	sim = sprintf("%d", 
								S[ coord ])
	
	err = sprintf("%d", 
								E[ coord ])
	
	list(sim=sim, err=err)
}


make.gen.matrix = function(break.idx, gen0, gen1, POS, S, E) {
	x = matrix(0, nrow = 3, ncol = 3, dimnames = list(c("0","1","2"), c("0","1","2")))
	
	state = get.gen.states(break.idx, gen0, gen1, S, E)
	
	for (k in 1:length(state$sim)) {
		#if (state$sim[k] != state$err[k]) {
		x[state$sim[k], state$err[k]] = x[state$sim[k], state$err[k]] + 1
		#}
	}
	#del = which(rowSums(x) == 0); if (length(del) > 0) x = x[-del, ]
	#del = which(colSums(x) == 0); if (length(del) > 0) x = x[, -del]
	x
}


plot.gen.matrix = function(mat, tag, name) {
	d = melt(mat)
	names(d) = c("Var1", "Var2", "x")
	
	d$Var1 = as.character(d$Var1)
	d$Var2 = as.character(d$Var2)
	
	s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]
	
	p = ggplot(data = d) + 
		geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
		geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
		scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme_bw() +
		theme(aspect.ratio = (nrow(mat))/((ncol(mat))*1),
					#axis.text.x = element_text(angle=-90),
					#legend.title = element_blank(),
					panel.border = element_blank(),
					#axis.ticks = element_blank(),
					panel.grid = element_blank()) +
		xlab("Breakpoint genotype, AFTER applying error profile") +
		ylab("Breakpoint genotype, BEFORE applying error profile") +
		ggtitle(toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag))))
	
	ggsave(filename = sprintf("_plot.break_generror.%s.%s.pdf", tag, name), plot = p, width = 10, height = 10, limitsize = FALSE)
}



get.gens.states = function(brk, g.pair0, g.pair1, S, E) {
	brk.pair0 = as.matrix(data.table(x = brk, y = g.pair0))
	brk.pair1 = as.matrix(data.table(x = brk, y = g.pair1))
	colnames(brk.pair0) = NULL
	colnames(brk.pair1) = NULL
	
	sim = sprintf("%d/%d", 
								S[ brk.pair0 ], 
								S[ brk.pair1 ])
	
	err = sprintf("%d/%d", 
								E[ brk.pair0 ], 
								E[ brk.pair1 ])
	
	list(sim=sim, err=err)
}


make.gens.matrix = function(break.idx, gen0, gen1, POS, S, E) {
	z = apply(expand.grid(c("0","1","2"), c("0","1","2")), 1, paste, collapse="/")
	x = matrix(0, nrow = length(z), ncol = length(z), dimnames = list(z, z))
	
	state = get.gens.states(break.idx, gen0, gen1, S, E)
	
	for (k in 1:length(state$sim)) {
		#if (state$sim[k] != state$err[k]) {
		x[state$sim[k], state$err[k]] = x[state$sim[k], state$err[k]] + 1
		#}
	}
	#del = which(rowSums(x) == 0); if (length(del) > 0) x = x[-del, ]
	#del = which(colSums(x) == 0); if (length(del) > 0) x = x[, -del]
	x
}


plot.gens.matrix = function(mat, tag, name, rm.nonexpect.cols = F) {
	if (rm.nonexpect.cols) {
		#u = expand.grid(c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"))
		#u = unlist(apply(u, 1, function(x) if (length(unique(x)) == 4) paste(x, collapse = ",") else NULL ))
		#u = sprintf("[%s]", u)
		u = c("0/2", "2/0")
		mat = mat[, which(colnames(mat) %in% u)]
	}
	
	d = melt(mat)
	names(d) = c("Var1", "Var2", "x")
	
	s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]
	
	p = ggplot(data = d) + 
		geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
		geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
		scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme_bw() +
		theme(aspect.ratio = (nrow(mat))/((ncol(mat))*1),
					#axis.text.x = element_text(angle=-90),
					#legend.title = element_blank(),
					panel.border = element_blank(),
					#axis.ticks = element_blank(),
					panel.grid = element_blank()) +
		xlab("Breakpoint genotypes, AFTER applying error profile") +
		ylab("Breakpoint genotypes, BEFORE applying error profile") +
		ggtitle(toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag))))
	
	ggsave(filename = sprintf("_plot.break_genotypes.%s.%s.pdf", tag, name), plot = p, width = 10, height = 10, limitsize = FALSE)
}




get.haps.states = function(brk, h.pair0.share, h.pair0.other, h.pair1.share, h.pair1.other, S, E) {
	brk.pair0.share = as.matrix(data.table(x = brk, y = h.pair0.share))
	brk.pair0.other = as.matrix(data.table(x = brk, y = h.pair0.other))
	brk.pair1.share = as.matrix(data.table(x = brk, y = h.pair1.share))
	brk.pair1.other = as.matrix(data.table(x = brk, y = h.pair1.other))
	colnames(brk.pair0.share) = NULL
	colnames(brk.pair0.other) = NULL
	colnames(brk.pair1.share) = NULL
	colnames(brk.pair1.other) = NULL
	
	sim = sprintf("%d%d/%d%d", 
								S[ brk.pair0.share ], S[ brk.pair0.other ], 
								S[ brk.pair1.share ], S[ brk.pair1.other ])
	
	err = sprintf("%d%d/%d%d", 
								E[ brk.pair0.share ], E[ brk.pair0.other ], 
								E[ brk.pair1.share ], E[ brk.pair1.other ])
	
	list(sim=sim, err=err)
}


make.haps.matrix = function(break.idx, hap0, hap1, POS, S, E) {
	z = apply(expand.grid(c("00","01","10","11"), c("00","01","10","11")), 1, paste, collapse="/")
	x = matrix(0, nrow = length(z), ncol = length(z), dimnames = list(z, z))
	
	h.pair0.share = hap0
	h.pair0.other = sapply(hap0, function(x) if (x %% 2 == 0) x - 1 else x + 1 )
	h.pair1.share = hap1
	h.pair1.other = sapply(hap1, function(x) if (x %% 2 == 0) x - 1 else x + 1 )
	
	state = get.haps.states(break.idx, h.pair0.share, h.pair0.other, h.pair1.share, h.pair1.other, S, E)
	
	for (k in 1:length(state$sim)) {
		#if (state$sim[k] != state$err[k]) {
		x[state$sim[k], state$err[k]] = x[state$sim[k], state$err[k]] + 1
		#}
	}
	#del = which(rowSums(x) == 0); if (length(del) > 0) x = x[-del, ]
	#del = which(colSums(x) == 0); if (length(del) > 0) x = x[, -del]
	x
}


plot.haps.matrix = function(mat, tag, name, rm.nonexpect.cols = F) {
	if (rm.nonexpect.cols) {
		#u = expand.grid(c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"))
		#u = unlist(apply(u, 1, function(x) if (length(unique(x)) == 4) paste(x, collapse = ",") else NULL ))
		#u = sprintf("[%s]", u)
		u = c("00/11", "01/10", "10/01", "11/00")
		mat = mat[, which(colnames(mat) %in% u)]
	}
	
	d = melt(mat)
	names(d) = c("Var1", "Var2", "x")
	
	s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]
	
	p = ggplot(data = d) + 
		geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
		geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
		scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme_bw() +
		theme(aspect.ratio = (nrow(mat))/((ncol(mat))*1),
					#axis.text.x = element_text(angle=-90),
					#legend.title = element_blank(),
					panel.border = element_blank(),
					#axis.ticks = element_blank(),
					panel.grid = element_blank()) +
		xlab("Breakpoint haplotypes, AFTER applying error profile") +
		ylab("Breakpoint haplotypes, BEFORE applying error profile") +
		ggtitle(toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag))))
	
	ggsave(filename = sprintf("_plot.break_haplotypes.%s.%s.pdf", tag, name), plot = p, width = 10, height = 10, limitsize = FALSE)
}



get.gamate.states = function(foc, brk, h.pair0.share, h.pair0.other, h.pair1.share, h.pair1.other, S, E) {
	brk.pair0.share = as.matrix(data.table(x = brk, y = h.pair0.share))
	brk.pair0.other = as.matrix(data.table(x = brk, y = h.pair0.other))
	brk.pair1.share = as.matrix(data.table(x = brk, y = h.pair1.share))
	brk.pair1.other = as.matrix(data.table(x = brk, y = h.pair1.other))
	colnames(brk.pair0.share) = NULL
	colnames(brk.pair0.other) = NULL
	colnames(brk.pair1.share) = NULL
	colnames(brk.pair1.other) = NULL
	
	foc.pair0.share = as.matrix(data.table(x = foc, y = h.pair0.share))
	foc.pair0.other = as.matrix(data.table(x = foc, y = h.pair0.other))
	foc.pair1.share = as.matrix(data.table(x = foc, y = h.pair1.share))
	foc.pair1.other = as.matrix(data.table(x = foc, y = h.pair1.other))
	colnames(foc.pair0.share) = NULL
	colnames(foc.pair0.other) = NULL
	colnames(foc.pair1.share) = NULL
	colnames(foc.pair1.other) = NULL
	
	sim = sprintf("[%d%d,%d%d,%d%d,%d%d]", 
								S[ foc.pair0.share ], S[ brk.pair0.share ], 
								S[ foc.pair0.other ], S[ brk.pair0.other ], 
								S[ foc.pair1.share ], S[ brk.pair1.share ], 
								S[ foc.pair1.other ], S[ brk.pair1.other ])
	
	err = sprintf("[%d%d,%d%d,%d%d,%d%d]", 
								E[ foc.pair0.share ], E[ brk.pair0.share ], 
								E[ foc.pair0.other ], E[ brk.pair0.other ], 
								E[ foc.pair1.share ], E[ brk.pair1.share ], 
								E[ foc.pair1.other ], E[ brk.pair1.other ])
	
	list(sim=sim, err=err)
}


make.gamate.matrix = function(focal.idx, break.idx, hap0, hap1, POS, S, E) {
	z = sprintf("[%s]", apply(expand.grid(c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11")), 1, paste, collapse=","))
	x = matrix(0, nrow = length(z), ncol = length(z), dimnames = list(z, z))
	
	h.pair0.share = hap0
	h.pair0.other = sapply(hap0, function(x) if (x %% 2 == 0) x - 1 else x + 1 )
	h.pair1.share = hap1
	h.pair1.other = sapply(hap1, function(x) if (x %% 2 == 0) x - 1 else x + 1 )
	
	state = get.gamate.states(focal.idx, break.idx, h.pair0.share, h.pair0.other, h.pair1.share, h.pair1.other, S, E)
	
	for (k in 1:length(state$sim)) {
		#if (state$sim[k] != state$err[k]) {
			x[state$sim[k], state$err[k]] = x[state$sim[k], state$err[k]] + 1
		#}
	}
	#del = which(rowSums(x) == 0); if (length(del) > 0) x = x[-del, ]
	#del = which(colSums(x) == 0); if (length(del) > 0) x = x[, -del]
	x
}


plot.gamete.matrix = function(mat, tag, name, rm.nonexpect.cols = F) {
	if (rm.nonexpect.cols) {
		#u = expand.grid(c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"), c("00","01","10","11"))
		#u = unlist(apply(u, 1, function(x) if (length(unique(x)) == 4) paste(x, collapse = ",") else NULL ))
		#u = sprintf("[%s]", u)
		u = c("[00,10,01,11]", "[00,11,01,10]", "[01,10,00,11]", "[01,11,00,10]", "[10,00,11,01]", "[10,01,11,00]", "[11,00,10,01]", "[11,01,10,00]")
		mat = mat[, which(colnames(mat) %in% u)]
	}
	
	d = melt(mat)
	names(d) = c("Var1", "Var2", "x")

	g = sort(unique(d$Var1))
	x = rbind(data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -1.5),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -1),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -0.5),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -0))
	
	g = sort(unique(d$Var2))
	y = rbind(data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -8),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -6),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -4),
						data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -2))

	s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]
	
	p = ggplot(data = d) + 
		geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
		geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
		geom_point(data=y, aes(x = gam, y = pos, colour=Gametes), fill = NA, shape=15, size=1) +
		geom_point(data=x, aes(x = pos, y = gam, colour=Gametes), fill = NA, shape=15, size=1) +
		scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
		scale_color_manual(values = c("0 | 0"="royalblue", "0 | 1"="limegreen", "1 | 0"="violetred1", "1 | 1"="orange")) +
		scale_x_discrete(expand = c(1,1)/100) +
		scale_y_discrete(expand = c(1,1)/100) +
		theme_bw() +
		theme(aspect.ratio = (nrow(mat)+5)/((ncol(mat)+10)*3),
					axis.text.x = element_text(angle=-90),
					#legend.title = element_blank(),
					panel.border = element_blank(),
					axis.ticks = element_blank(),
					panel.grid = element_blank(),
					axis.text = element_blank()) +
		xlab("Breakpoint gametic states, AFTER applying error profile") +
		ylab("Breakpoint gametic states, BEFORE applying error profile") +
		ggtitle(toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag))))
	
	ggsave(filename = sprintf("_plot.fourgamates.%s.%s.pdf", tag, name), plot = p, width = ncol(mat)/20, height = nrow(mat)/20, limitsize = FALSE)
}


err = split(data, data$origin)

sim = err$simulated
err$simulated = NULL

for (tag in names(err)) {
	a = which(err[[tag]]$g.end.lhs | err[[tag]]$g.end.rhs)
	b = which(err[[tag]]$h.end.lhs | err[[tag]]$h.end.rhs)
	c = which(err[[tag]]$p.end.lhs | err[[tag]]$p.end.rhs)
	del = unique(c(a, b, c))
	if (length(del) > 0) {
		err[[tag]] = err[[tag]][-del, ]
	}
}

SG = load.bigmatrix(sprintf("%s.G", prefix))
SH = load.bigmatrix(sprintf("%s.H", prefix))
SP = load.bigmatrix(sprintf("%s.P", prefix))

root = getwd()

err.gen = list()
err.gens = list()
err.haps = list()
err.haps.gam = list()
err.phas = list()
err.phas.gam = list()

for (tag in names(err)) {
	cat("Error detection in:", tag, "\n")

	sub = err[[tag]]
	
	setwd(tag)
	load(sprintf("%s.%s.RData", prefix, tag))
	EG = load.bigmatrix(sprintf("%s.%s.G", prefix, tag))
	EH = load.bigmatrix(sprintf("%s.%s.H", prefix, tag))
	EP = load.bigmatrix(sprintf("%s.%s.P", prefix, tag))
	setwd(root)
	
	
	
	# genotypes
	cat(" genotypes\n")
	brk = c(match(sub$g.lhs, POS) - 1,
					match(sub$g.rhs, POS) + 1)
	g0 = c(sub$g0, sub$g0)
	g1 = c(sub$g1, sub$g1)
	
	err.gen[[tag]] = make.gen.matrix(brk, g0, g1, POS, SG, EG)
	#mat = mat / sum(mat)
	#plot.gen.matrix(mat, tag, "genotypes")
	
	err.gens[[tag]] = make.gens.matrix(brk, g0, g1, POS, SG, EG)
	#mat = mat / sum(mat)
	#plot.gens.matrix(mat, tag, "genotypes")
	

	foc = match(sub$position, POS)
	foc = c(foc, foc)
	
	# haplotypes
	cat(" haplotypes\n")
	brk = c(match(sub$h.lhs, POS) - 1,
					match(sub$h.rhs, POS) + 1)
	h0 = c(sub$h0, sub$h0)
	h1 = c(sub$h1, sub$h1)
	
	err.haps[[tag]] = make.haps.matrix(brk, h0, h1, POS, SH, EH)
	#mat = mat / sum(mat)
	#plot.haps.matrix(mat, tag, "haplotypes")
	
	err.haps.gam[[tag]] = make.gamate.matrix(foc, brk, h0, h1, POS, SH, EH)
	#mat = mat / sum(mat)
	#plot.gamete.matrix(mat, tag, "haplotypes", rm.nonexpect.cols = T)
	
	# phased haplotypes
	cat(" phased haplotypes\n")
	brk = c(match(sub$p.lhs, POS) - 1,
					match(sub$p.rhs, POS) + 1)
	p0 = c(sub$p0, sub$p0)
	p1 = c(sub$p1, sub$p1)
	
	err.phas[[tag]] = make.haps.matrix(brk, p0, p1, POS, SP, EP)
	#mat = mat / sum(mat)
	#plot.haps.matrix(mat, tag, "phased_haplotypes")
	
	err.phas.gam[[tag]] = make.gamate.matrix(foc, brk, p0, p1, POS, SP, EP)
	#mat = mat / sum(mat)
	#plot.gamete.matrix(mat, tag, "phased_haplotypes", rm.nonexpect.cols = T)

}



save(err.gen,
		 err.gens,
		 err.haps,
		 err.haps.gam,
		 err.phas,
		 err.phas.gam,
		 file = "_plotdata.break.RData")

#load("_plotdata.break.RData")


# gen

for (tag in names(err.gen)) {
	err.gen[[tag]] = melt(err.gen[[tag]])
	err.gen[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.gen)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=1.5) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 16),
				#axis.text.x = element_text(angle=-90),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original breakpoint genotype") +
	xlab("Breakpoint genotype, after applying error profile") +
	ggtitle("Genotype")

ggsave(filename = "_plot.break_gen.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)




# gens

for (tag in names(err.gens)) {
	err.gens[[tag]] = melt(err.gens[[tag]])
	err.gens[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.gens)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=1) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 12),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original breakpoint genotype pair") +
	xlab("Breakpoint genotype pair, after applying error profile") +
	ggtitle("Genotype pair")

ggsave(filename = "_plot.break_gen_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)



# haps

for (tag in names(err.haps)) {
	err.haps[[tag]] = melt(err.haps[[tag]])
	err.haps[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.haps)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 10),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original breakpoint haplotype pair") +
	xlab("Breakpoint haplotype pair, after applying error profile") +
	ggtitle("Haplotype pair")

ggsave(filename = "_plot.break_hap_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)



# phas

for (tag in names(err.phas)) {
	err.phas[[tag]] = melt(err.phas[[tag]])
	err.phas[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.phas)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 10),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original breakpoint phased haplotype pair") +
	xlab("Breakpoint phased haplotype pair, after applying error profile") +
	ggtitle("Phased haplotype pair")

ggsave(filename = "_plot.break_phap_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)



# haps.gam

u = c("[00,10,01,11]", "[00,11,01,10]", "[01,10,00,11]", "[01,11,00,10]", "[10,00,11,01]", "[10,01,11,00]", "[11,00,10,01]", "[11,01,10,00]")

for (tag in names(err.haps.gam)) {
	err.haps.gam[[tag]] = err.haps.gam[[tag]][, which(colnames(err.haps.gam[[tag]]) %in% u)]
	err.haps.gam[[tag]] = melt(err.haps.gam[[tag]])
	err.haps.gam[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.haps.gam)
names(d) = c("Var1", "Var2", "x", "profile")

g = sort(unique(d$Var1))
x = rbind(data.table(gam = g, Gametes = "---", pos = -1.75),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -1.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -1),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -0.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -0))

g = sort(unique(d$Var2))
y = rbind(data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -6),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -4.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -3),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -1.5))

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	geom_point(data=y, aes(x = gam, y = pos, colour=Gametes), fill = NA, shape=15, size=1) +
	geom_point(data=x, aes(x = pos, y = gam, colour=Gametes), fill = NA, shape=15, size=1) +
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_color_manual(values = c("0 | 0"="royalblue", "0 | 1"="limegreen", "1 | 0"="violetred1", "1 | 1"="orange", "---"="white")) +
	scale_x_discrete(expand = c(1,1)/200) +
	scale_y_discrete(expand = c(1,1)/200) +
	theme_classic() +
	theme(aspect.ratio = (256+2)/(8*5),
				axis.text = element_blank(),
				#axis.text.x = element_text(angle=-90),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	xlab("Breakpoint gametic state, after applying error profile") +
	ylab("Original breakpoint gametic state") +
	ggtitle("Gametes of haplotypes")

ggsave(filename = "_plot.gamates_hap.pdf", plot = p, width = 15, height = 15, limitsize = FALSE)




# phas.gam

u = c("[00,10,01,11]", "[00,11,01,10]", "[01,10,00,11]", "[01,11,00,10]", "[10,00,11,01]", "[10,01,11,00]", "[11,00,10,01]", "[11,01,10,00]")

for (tag in names(err.phas.gam)) {
	err.phas.gam[[tag]] = err.phas.gam[[tag]][, which(colnames(err.phas.gam[[tag]]) %in% u)]
	err.phas.gam[[tag]] = melt(err.phas.gam[[tag]])
	err.phas.gam[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, err.phas.gam)
names(d) = c("Var1", "Var2", "x", "profile")

g = sort(unique(d$Var1))
x = rbind(data.table(gam = g, Gametes = "---", pos = -1.75),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -1.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -1),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -0.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -0))

g = sort(unique(d$Var2))
y = rbind(data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\1 | \\2", g), pos = -6),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\3 | \\4", g), pos = -4.5),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\5 | \\6", g), pos = -3),
					data.table(gam = g, Gametes = sub("^\\[([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1]),([0-1])([0-1])\\]$", "\\7 | \\8", g), pos = -1.5))

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	geom_point(data=y, aes(x = gam, y = pos, colour=Gametes), fill = NA, shape=15, size=1) +
	geom_point(data=x, aes(x = pos, y = gam, colour=Gametes), fill = NA, shape=15, size=1) +
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_color_manual(values = c("0 | 0"="royalblue", "0 | 1"="limegreen", "1 | 0"="violetred1", "1 | 1"="orange", "---"="white")) +
	scale_x_discrete(expand = c(1,1)/200) +
	scale_y_discrete(expand = c(1,1)/200) +
	theme_classic() +
	theme(aspect.ratio = (256+2)/(8*5),
				axis.text = element_blank(),
				#axis.text.x = element_text(angle=-90),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	xlab("Breakpoint gametic state, after applying error profile") +
	ylab("Original breakpoint gametic state") +
	ggtitle("Gametes of phased haplotypes")

ggsave(filename = "_plot.gamates_phap.pdf", plot = p, width = 15, height = 15, limitsize = FALSE)






##################################################
##################################################
##################################################

### inspect missed RVs

err = split(data, data$origin)

for (tag in names(err)) {
	a = which(err[[tag]]$g.end.lhs | err[[tag]]$g.end.rhs)
	b = which(err[[tag]]$h.end.lhs | err[[tag]]$h.end.rhs)
	c = which(err[[tag]]$p.end.lhs | err[[tag]]$p.end.rhs)
	del = unique(c(a, b, c))
	if (length(del) > 0) {
		err[[tag]] = err[[tag]][-del, ]
	}
}

sim = err$simulated
err$simulated = NULL


SG = load.bigmatrix(sprintf("%s.G", prefix))
SH = load.bigmatrix(sprintf("%s.H", prefix))
SP = load.bigmatrix(sprintf("%s.P", prefix))

root = getwd()

missed.gen = list()
missed.gens = list()
missed.haps = list()
missed.haps.gam = list()
missed.phas = list()
missed.phas.gam = list()

for (tag in names(err)) {
	cat("Error detection in:", tag, "\n")
	
	setwd(tag)
	load(sprintf("%s.%s.RData", prefix, tag))
	EG = load.bigmatrix(sprintf("%s.%s.G", prefix, tag))
	EH = load.bigmatrix(sprintf("%s.%s.H", prefix, tag))
	EP = load.bigmatrix(sprintf("%s.%s.P", prefix, tag))
	setwd(root)
	
	
	foc = match(sim$position, POS)

	
	# genotypes
	cat(" genotypes\n")
	missed.gen[[tag]] = make.gen.matrix(foc, sim$g0, sim$g1, POS, SG, EG)
	#mat = mat / sum(mat)
	#plot.gen.matrix(mat, tag, "genotypes")
	
	missed.gens[[tag]] = make.gens.matrix(foc, sim$g0, sim$g1,POS, SG, EG)
	#mat = mat / sum(mat)
	#plot.gens.matrix(mat, tag, "genotypes")
	
	

	
	# haplotypes
	cat(" haplotypes\n")
	missed.haps[[tag]] = make.haps.matrix(foc, sim$h0, sim$h1, POS, SH, EH)
	#mat = mat / sum(mat)
	#plot.haps.matrix(mat, tag, "haplotypes")

	
	# phased haplotypes
	cat(" phased haplotypes\n")
	missed.phas[[tag]] = make.haps.matrix(foc, sim$p0, sim$p1, POS, SP, EP)
	#mat = mat / sum(mat)
	#plot.haps.matrix(mat, tag, "phased_haplotypes")

}




save(missed.gen,
		 missed.gens,
		 missed.haps,
		 missed.phas,
		 file = "_plotdata.missing.RData")

#load("_plotdata.missing.RData")


# gen

for (tag in names(missed.gen)) {
	missed.gen[[tag]] = melt(missed.gen[[tag]])
	missed.gen[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, missed.gen)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=1.5) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 16),
				#axis.text.x = element_text(angle=-90),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original focal genotype") +
	xlab("Focal genotype, after applying error profile") +
	ggtitle("Genotype")

ggsave(filename = "_plot.missing_gen.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)




# gens

for (tag in names(missed.gens)) {
	missed.gens[[tag]] = melt(missed.gens[[tag]])
	missed.gens[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, missed.gens)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=1) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 12),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original focal genotype pair") +
	xlab("Focal genotype pair, after applying error profile") +
	ggtitle("Genotype pair")

ggsave(filename = "_plot.missing_gen_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)



# haps

for (tag in names(missed.haps)) {
	missed.haps[[tag]] = melt(missed.haps[[tag]])
	missed.haps[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, missed.haps)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 10),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original focal haplotype pair") +
	xlab("Focal haplotype pair, after applying error profile") +
	ggtitle("Haplotype pair")

ggsave(filename = "_plot.missing_hap_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)



# phas

for (tag in names(missed.phas)) {
	missed.phas[[tag]] = melt(missed.phas[[tag]])
	missed.phas[[tag]]$profile = toupper(gsub("_", " ", sub("^[a-z]+_(.+)$", "\\1", tag)))
}
d = Reduce(rbind, missed.phas)
names(d) = c("Var1", "Var2", "x", "profile")

d$Var1 = as.character(d$Var1)
d$Var2 = as.character(d$Var2)

s = d[which(apply(d, 1, function(x) (x[1] == x[2] && as.numeric(x[3]) != 0) )), ]

p = ggplot(data = d) + 
	facet_wrap(~profile, nrow=1) +
	geom_tile(aes(x = Var2, y = Var1, fill = log10(x)), colour = "grey50") + 
	geom_tile(data=s, aes(x = Var2, y = Var1), colour = "green2", fill="black", alpha=0, size=0.75) + 
	scale_fill_gradient(low = "white", high = "darkred", na.value = "grey85") +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.text = element_text(family = "Courier", size = 10),
				axis.text.x = element_text(angle=-90, vjust = 0.5),
				#legend.title = element_blank(),
				panel.border = element_blank(),
				#axis.ticks = element_blank(),
				panel.grid = element_blank()) +
	ylab("Original focal phased haplotype pair") +
	xlab("Focal phased haplotype pair, after applying error profile") +
	ggtitle("Phased haplotype pair")

ggsave(filename = "_plot.missing_phap_pair.pdf", plot = p, width = 14, height = 4, limitsize = FALSE)














