#
# plot error profile
#

library(ggplot2)
library(grid)
library(data.table)


args = commandArgs(T)

pro.file = args[1]
par.rm.nontyped = as.logical(args[2])


named = sub("^profile\\.(.+)\\.RData$", "\\1", pro.file)
typed = as.character(par.rm.nontyped)


cat("Profile name:", named, "\n")
cat("Include typed only:", typed, "\n")


cat("Loading ... ")
load(pro.file)
cat("OK\n")



cat("Calculating allele frequency spectrum ... ")
afs = by(p, p$maf, function(x) {
	data.frame(maf = x$maf[1], count = nrow(x))
})
afs = rbindlist(afs)
cat("OK\n")




cat("Filtering ... ")

filter.profile.data = function(p, rm.nontyped = F) {
	del = which(is.na(p$true.gt))
	if (length(del) > 0) {
		p = p[-del, ]
	}
	
	del = which(is.na(p$call.gt))
	if (length(del) > 0) {
		p = p[-del, ]
	}
	
	p$error = (p$true.gt != p$call.gt)
	
	del = which(p$maf == 0 & !p$error)
	if (length(del) > 0) {
		p = p[-del, ]
	}
	
	del = which(! p$true.is.confident)
	if (length(del) > 0) {
		p = p[-del, ]
	}
	
	if (rm.nontyped) {
		# 		del = which(! p$true.is.typed)
		# 		if (length(del) > 0) {
		# 			p = p[-del, ]
		# 		}
		del = which(! p$call.is.typed)
		if (length(del) > 0) {
			p = p[-del, ]
		}
	}
	
	p
}

p = filter.profile.data(p, par.rm.nontyped)

cat("OK\n")




cat("Making alt allele the minor allele ... ")
i = which(! p$alt.is.minor)
if (length(i) > 0) {
	j = which(p$true.gt[i] == 0); p$true.gt[i][j] = abs(p$true.gt[i][j] - 2)
	j = which(p$call.gt[i] == 0); p$call.gt[i][j] = abs(p$call.gt[i][j] - 2)
}
cat("OK\n")




cat("Re-calculating allele frequency spectrum ... ")
afs.post = by(p, p$maf, function(x) {
	data.frame(maf = x$maf[1], count = nrow(x))
})
afs.post = rbindlist(afs.post)
cat("OK\n")


d.afs = rbind(cbind(afs, label = "Reference"),
							cbind(afs.post, label = "Available"))


gg = ggplot(data = d.afs[d.afs$maf != 0, ]) + 
	geom_point(aes(x=maf, y=count, colour = label)) + 
	scale_colour_manual(values = c("Reference" = "black", "Available" = "peru")) +
	scale_x_log10(expand = c(0.01, 0.01), breaks = c(0.05, 0.5, 5, 50) / 100, labels = c(0.05, 0.5, 5, 50)) + 
	scale_y_log10(expand = c(0.01, 0.01), breaks = 10^(1:8), labels = sprintf("%d", 10^(1:8))) +
	coord_cartesian(xlim = c(0.0001, 0.5)) +
	theme_classic() +
	theme(panel.grid.major = element_line(colour="grey90"),
				legend.title = element_blank(),
				aspect.ratio = 1) +
	xlab("Minor allele frequency (%), log-scale") +
	ylab("Count, log-scale")

ggsave(sprintf("_plot.afs.typed_%s.pdf", named, typed), gg, width = 10, height = 10)



cat("Labelling ... ")

p$label.gt = NA

p$label.gt[ which(p$true.gt == 0 & p$call.gt == 0) ] = "Truth 0, typed 0" # correct
p$label.gt[ which(p$true.gt == 0 & p$call.gt == 1) ] = "Truth 0, typed 1"
p$label.gt[ which(p$true.gt == 0 & p$call.gt == 2) ] = "Truth 0, typed 2"

p$label.gt[ which(p$true.gt == 1 & p$call.gt == 0) ] = "Truth 1, typed 0"
p$label.gt[ which(p$true.gt == 1 & p$call.gt == 1) ] = "Truth 1, typed 1" # correct
p$label.gt[ which(p$true.gt == 1 & p$call.gt == 2) ] = "Truth 1, typed 2"

p$label.gt[ which(p$true.gt == 2 & p$call.gt == 0) ] = "Truth 2, typed 0"
p$label.gt[ which(p$true.gt == 2 & p$call.gt == 1) ] = "Truth 2, typed 1"
p$label.gt[ which(p$true.gt == 2 & p$call.gt == 2) ] = "Truth 2, typed 2" # correct

if (any(! p$call.is.typed)) p$label.gt[ which(! p$call.is.typed) ] = "Not typed"
p$label.gt = factor(p$label.gt)


p$label.true.gt = NA

p$label.true.gt[ which(p$true.gt == 0) ] = "Truth genotype = 0"
p$label.true.gt[ which(p$true.gt == 1) ] = "Truth genotype = 1"
p$label.true.gt[ which(p$true.gt == 2) ] = "Truth genotype = 2"

p$label.true.gt = factor(p$label.true.gt)


p$label.call.gt = NA

p$label.call.gt[ which(p$call.gt == 0) ] = "Typed genotype = 0"
p$label.call.gt[ which(p$call.gt == 1) ] = "Typed genotype = 1"
p$label.call.gt[ which(p$call.gt == 2) ] = "Typed genotype = 2"

if (any(! p$call.is.typed)) p$label.call.gt[ which(! p$call.is.typed) ] = "Not typed"
p$label.call.gt = factor(p$label.call.gt)


p$type = NA

p$type[ which(p$true.gt == 0 & p$call.gt == 0) ] = "00" # correct
p$type[ which(p$true.gt == 0 & p$call.gt == 1) ] = "01"
p$type[ which(p$true.gt == 0 & p$call.gt == 2) ] = "02"

p$type[ which(p$true.gt == 1 & p$call.gt == 0) ] = "10"
p$type[ which(p$true.gt == 1 & p$call.gt == 1) ] = "11" # correct
p$type[ which(p$true.gt == 1 & p$call.gt == 2) ] = "12"

p$type[ which(p$true.gt == 2 & p$call.gt == 0) ] = "20"
p$type[ which(p$true.gt == 2 & p$call.gt == 1) ] = "21"
p$type[ which(p$true.gt == 2 & p$call.gt == 2) ] = "22" # correct

p$type = factor(p$type)




q = p[which(p$chr==1), ]

q$state.true = NA
q$state.true[ which(q$true.gt != 1) ] = "hom"
q$state.true[ which(q$true.gt == 1) ] = "het"

q$state.call = NA
q$state.call[ which(q$call.gt != 1) ] = "hom"
q$state.call[ which(q$call.gt == 1) ] = "het"

q$gt0 = q$ref.frq^2
q$gt1 = 2 * q$ref.frq * q$alt.frq
q$gt2 = q$alt.frq^2

q$gt.true = NA
i = which(q$true.gt == 0); q$gt.true[i] = q$gt0[i]
i = which(q$true.gt == 1); q$gt.true[i] = q$gt1[i]
i = which(q$true.gt == 2); q$gt.true[i] = q$gt2[i]

q$gt.call = NA
i = which(q$call.gt == 0); q$gt.call[i] = q$gt0[i]
i = which(q$call.gt == 1); q$gt.call[i] = q$gt1[i]
i = which(q$call.gt == 2); q$gt.call[i] = q$gt2[i]


i = sort(sample(nrow(q), 100000))
ggplot(data = q[which(q$error), ]) + geom_point(aes(x=(gt.true), y=(gt.call), colour=as.character(abs(call.gt-true.gt)))) + theme(aspect.ratio=1) + scale_x_log10() + scale_y_log10()
ggplot(data = q[which(q$error), ]) + geom_point(aes(x=gt.true + (gt.true - (1-gt.call)), y=gt.call + (gt.call - (1-gt.true)), colour=interaction(state.true, state.call))) + theme(aspect.ratio=1)
ggplot(data = q[which(q$error), ]) + geom_point(aes(x=gt.true, y=(gt.true - gt.call), colour=interaction(state.true, state.call))) + theme(aspect.ratio=1)
ggplot(data = q[which(q$error), ]) + facet_grid(label.call.gt~label.true.gt) + geom_point(aes(x=(gt.true), y=(gt.call))) + theme(aspect.ratio=1)



load(sprintf("_stat.prop_matrix.%s.typed_%s.RData", named, typed))
q$emp = NA

q$emp[ which(q$true.gt == 0 & q$call.gt == 0) ] = m[1, 1]
q$emp[ which(q$true.gt == 0 & q$call.gt == 1) ] = m[1, 2]
q$emp[ which(q$true.gt == 0 & q$call.gt == 2) ] = m[1, 3]

q$emp[ which(q$true.gt == 1 & q$call.gt == 0) ] = m[2, 1]
q$emp[ which(q$true.gt == 1 & q$call.gt == 1) ] = m[2, 2]
q$emp[ which(q$true.gt == 1 & q$call.gt == 2) ] = m[2, 3]

q$emp[ which(q$true.gt == 2 & q$call.gt == 0) ] = m[3, 1]
q$emp[ which(q$true.gt == 2 & q$call.gt == 1) ] = m[3, 2]
q$emp[ which(q$true.gt == 2 & q$call.gt == 2) ] = m[3, 3]


ggplot(data = q[which(q$error), ]) + geom_point(aes(x= ref.frq, y= alt.frq * emp, colour = type)) + theme(aspect.ratio=1)



predict.error = function(true, ref, alt) {
	sample(c(0, 1, 2), 1, TRUE, c(ref^2, 2*ref*alt, alt^2))
}

q$pred.gt = apply(q, 1, function(x) predict.error(x$true.gt, x$ref.frq, x$alt.frq))



d = rbindlist(by(p, p$maf, function(x) {
	z = rbindlist(by(x, (x$label.gt), function(y) {
		if (nrow(y) == 0) return(NULL)
		z = y[1, ]
		z$count.type = nrow(y)
		z
	}))
	z$count.maf = nrow(x)
	z
}))

ggplot(data=d) + facet_wrap(~label.true.gt, nrow=1) + geom_bar(aes(x=maf, y=count.type, fill=label.call.gt), stat="identity", position="fill")

ggplot(data=d) + geom_line(aes(x=maf, y=count.type, colour=label.gt)) + scale_y_log10()

ggplot(data = p) + geom_freqpoly(aes(x=maf, colour = abs(call.gt-true.gt)))



# MAF bin breaks
brk = seq(0, 0.5, length.out = 101) # c(0, exp(seq(log(0.0005), log(0.5), length.out = 100))) # 
p$bin = cut(p$maf, breaks = brk, include.lowest = T)

cat("OK\n")



cat("Preparing plot data ... ")

make.data.prop.together = function(p) {
	splt = split(p, p$bin)
	splt = lapply(splt, function(x) {
		if (nrow(x) == 0) return(NULL)
		sub = as.data.table(x)
		sub = split(sub, sub$type)
		sub = lapply(sub, function(y, n) {
			if (nrow(y) == 0) return(NULL)
			z = y[1, ]
			z$count = nrow(y)
			z$percent = nrow(y) / n * 100
			z
		}, nrow(x))
		sub = rbindlist(sub)
		if (nrow(sub) == 0) return(NULL)
		sub$nbin = nrow(x)
		sub
	})
	as.data.frame(rbindlist(splt))
}


make.data.prop.bytruegt = function(p) {
	splt = split(p, list(p$bin, p$true.gt))
	splt = lapply(splt, function(x) {
		if (nrow(x) == 0) return(NULL)
		sub = as.data.table(x)
		sub = split(sub, sub$type)
		sub = lapply(sub, function(y, n) {
			if (nrow(y) == 0) return(NULL)
			z = y[1, ]
			z$count = nrow(y)
			z$percent = nrow(y) / n * 100
			z
		}, nrow(x))
		sub = rbindlist(sub)
		if (nrow(sub) == 0) return(NULL)
		sub$nbin = nrow(x)
		sub
	})
	as.data.frame(rbindlist(splt))
}


make.data.pos.together = function(p) {
	rng = range(p$pos)
	brk = seq(rng[1], rng[2], length.out = 501)
	inv =  cut(p$pos, breaks = brk, include.lowest = T)
	p$inv = inv
	
	splt = split(p, list(p$chr, p$inv))
	splt = lapply(splt, function(x) {
		if (nrow(x) == 0) return(NULL)
		z = split(x, x$type)
		z = lapply(z, function(y, n) {
			if (nrow(y) == 0) return(NULL)
			data.table(y[1, ],
								 stat.count = nrow(y),
								 stat.percent = nrow(y) / n * 100)
		}, nrow(x))
		z = rbindlist(z)
		z$inv.size = nrow(x)
		z
	})
	d = as.data.frame(rbindlist(splt))
	d$inv.beg = as.numeric(sub("^.{1}(.+),(.+).{1}$", "\\1", as.character(d$inv)))/1e06
	d$inv.end = as.numeric(sub("^.{1}(.+),(.+).{1}$", "\\2", as.character(d$inv)))/1e06
	d$inv.mid = d$inv.beg + ((d$inv.end - d$inv.beg) / 2)
	d
}


make.data.pos.bytruegt = function(p) {
	rng = range(p$pos)
	brk = seq(rng[1], rng[2], length.out = 501)
	inv =  cut(p$pos, breaks = brk, include.lowest = T)
	p$inv = inv
	
	splt = split(p, list(p$chr, p$inv, p$true.gt))
	splt = lapply(splt, function(x) {
		if (nrow(x) == 0) return(NULL)
		z = split(x, x$type)
		z = lapply(z, function(y, n) {
			if (nrow(y) == 0) return(NULL)
			data.table(y[1, ],
								 stat.count = nrow(y),
								 stat.percent = nrow(y) / n * 100)
		}, nrow(x))
		z = rbindlist(z)
		z$inv.size = if (nrow(z) == 0) 0 else nrow(x)
		z
	})
	d = as.data.frame(rbindlist(splt))
	d$inv.beg = as.numeric(sub("^.{1}(.+),(.+).{1}$", "\\1", as.character(d$inv)))/1e06
	d$inv.end = as.numeric(sub("^.{1}(.+),(.+).{1}$", "\\2", as.character(d$inv)))/1e06
	d$inv.mid = d$inv.beg + ((d$inv.end - d$inv.beg) / 2)
	d
}


d.prop.together = make.data.prop.together(p)
d.prop.bytruegt = make.data.prop.bytruegt(p)
d.pos.together = make.data.pos.together(p)
d.pos.bytruegt = make.data.pos.bytruegt(p)


save(d.afs,
		 d.prop.together,
		 d.prop.bytruegt,
		 d.pos.together,
		 d.pos.bytruegt,
		 file = sprintf("_data.%s.typed_%s.RData", pro.file, as.character(par.rm.nontyped)))


cat("OK\n")


#
# plot 
#

cat("Plotting ... ")

colour.gt = c("Truth 0, called 0" = "peachpuff",
							"Truth 0, called 1" = "orange",
							"Truth 0, called 2" = "orangered",
							"Truth 1, called 0" = "green2",
							"Truth 1, called 1" = "darkseagreen1",
							"Truth 1, called 2" = "seagreen",
							"Truth 2, called 0" = "dodgerblue4",
							"Truth 2, called 1" = "deepskyblue",
							"Truth 2, called 2" = "lightblue1")

colour.true.gt = c("Truth genotype = 0" = "orange",
									 "Truth genotype = 1" = "green2",
									 "Truth genotype = 2" = "deepskyblue")

linety.true.gt = c("Truth genotype = 0" = "solid",
									 "Truth genotype = 1" = "dashed",
									 "Truth genotype = 2" = "dotted")

colour.call.gt = c("Called genotype = 0" = "sienna1",
									 "Called genotype = 1" = "olivedrab",
									 "Called genotype = 2" = "mediumpurple3")

linety.call.gt = c("Called genotype = 0" = "solid",
									 "Called genotype = 1" = "dashed",
									 "Called genotype = 2" = "dotted")



ggplot(data = p) + geom_bar(aes(x=as.numeric(bin)-0.5, fill=label.true.gt), position = "fill") 

ggplot(data = p) + geom_bar(aes(x=as.numeric(bin)-0.5, fill=label.call.gt), position = "fill") 



d = by(p, p$maf, function(x) {
	z = x[1, ]
	z$N = nrow(x)
	z
})
d = rbindlist(d)
ggplot(data = d) + geom_point(aes(x=maf, y=N)) + scale_x_log10(breaks=seq(5, 50, by=5)/100) + scale_y_log10()


d = by(p, p$bin, function(x) {
	z = x[1, ]
	z$bin.error.prop = length(which(x$error)) / nrow(x)
	z
})
d = rbindlist(d)

ggplot(data = d) + geom_line(aes(x=as.numeric(bin)-0.5, y = bin.error.prop))



plot.prop.together = function(d, ymax = 100, yticks = 21) {
	ggplot(data = d) + 
		geom_bar(aes(x=as.numeric(bin)-0.5, y=percent, fill=label.gt), stat="identity") + 
		scale_fill_manual(values = colour.gt) +
		scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
		scale_y_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, ymax, length.out = yticks)) +
		coord_cartesian(ylim = c(0, ymax)) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					legend.title = element_blank()) +
		xlab("MAF bin (%)") +
		ylab("Relative Proportion (%)")
}


plot.size.together = function(d) {
	n = data.frame(bin = d$bin, nbin = d$nbin / 1e03)
	n = unique(n)
	if (any(n$nbin < 1)) {
		n$nbin[which(n$nbin < 1)] = 1
	}
	
	ggplot(data = n) + 
		geom_bar(aes(x=as.numeric(bin)-0.5, y=nbin), alpha=0.5, stat="identity") + 
		scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
		scale_y_log10(expand = c(0.0025, 0.0025), breaks = 10^(0:8) ) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					legend.title = element_blank()) +
		xlab("MAF bin (%)") +
		ylab("Number of markers (thousands), log-scale")
}


plot.prop.bytruegt = function(d, ymax = 100, yticks = 21) {
	ggplot(data = d) + 
		facet_wrap(~label.true.gt) +
		geom_bar(aes(x=as.numeric(bin)-0.5, y=percent, fill=label.gt), stat="identity") + 
		scale_fill_manual(values = colour.gt) +
		scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
		scale_y_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, ymax, length.out = yticks)) +
		coord_cartesian(ylim = c(0, ymax)) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					legend.title = element_blank(),
					strip.background = element_rect(fill = NA, colour = NA),
					panel.border = element_rect(fill = NA, colour = "black"),
					panel.margin = unit(0.25, units = "inches")) +
		xlab("MAF bin (%)") +
		ylab("Relative Proportion (%)")
}


plot.prop.bytruegt.lines = function(d, ymax = 100, yticks = 21) {
	ggplot(data = d) + 
		facet_wrap(~label.true.gt) +
		geom_line(aes(x=as.numeric(bin)-0.5, y=percent, colour=label.gt), size=1) + 
		scale_colour_manual(values = colour.gt) +
		scale_x_continuous(expand = c(0.005, 0.005), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
		scale_y_continuous(expand = c(0.005, 0.005), breaks = seq(0, ymax, length.out = yticks)) +
		coord_cartesian(ylim = c(0, ymax)) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					legend.title = element_blank(),
					strip.background = element_rect(fill = NA, colour = NA),
					panel.border = element_rect(fill = NA, colour = "black"),
					panel.margin = unit(0.25, units = "inches")) +
		xlab("MAF bin (%)") +
		ylab("Relative Proportion (%)")
}


plot.size.bytruegt = function(d) {
	n = data.frame(bin = d$bin, nbin = d$nbin, label.true.gt = as.character(d$label.true.gt))
	n = unique(n)
	if (any(n$nbin < 1)) {
		n$nbin[which(n$nbin < 1)] = 1
	}
	
	ggplot(data = n) + 
		facet_wrap(~label.true.gt) +
		geom_bar(aes(x=as.numeric(bin)-0.5, y=(nbin+1)), alpha=0.5, stat="identity", position = "identity") + 
		scale_fill_manual(values = colour.gt) +
		scale_x_continuous(expand = c(0.0025, 0.0025), breaks = seq(0, 100, length.out = 11), labels = seq(0, 50, length.out = 11)) +
		scale_y_log10(expand = c(0.0025, 0.0025), breaks = (10^(0:10)) ) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					legend.title = element_blank(),
					strip.background = element_rect(fill = NA, colour = NA),
					panel.border = element_rect(fill = NA, colour = "black"),
					panel.margin = unit(0.25, units = "inches")) +
		xlab("MAF bin (%)") +
		ylab("Number of markers (thousands), log-scale")
}


plot.pos.together = function(d) {
	ggplot(data = d) + 
		facet_wrap(~chr, ncol = 2) +
		geom_rect(data=d, aes(xmin=inv.beg, xmax=inv.end, ymin=-Inf, ymax=Inf), colour="grey90", alpha=0.05) +
		geom_line(aes(x=inv.mid, y=stat.percent, colour=label.call.gt)) + 
		#geom_bar(aes(x=inv, y=stat.percent, fill=label.gt), stat="identity", position="fill") + 
		scale_colour_manual(values = colour.call.gt) +
		scale_x_continuous(expand = c(0.005, 0.005)) +
		#coord_cartesian(ylim = c(-7.5, 102.5)) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					strip.text = element_text(face = "bold"),
					legend.title = element_blank(),
					strip.background = element_rect(fill = NA, colour = NA),
					panel.border = element_rect(fill = NA, colour = "black"),
					panel.margin = unit(0.25, units = "inches"),
					panel.margin.y = unit(0.1, units = "inches"),
					legend.position = "bottom") +
		xlab("Chromosome position (Mb)") +
		ylab("Relative error proportion (%)")
}


plot.pos.bytruegt = function(d) {
	ggplot(data = d) + 
		facet_grid(chr~label.true.gt) +
		geom_rect(data=d, aes(xmin=inv.beg, xmax=inv.end, ymin=-Inf, ymax=Inf), colour="grey90", alpha=0.05) +
		geom_line(aes(x=inv.mid, y=stat.percent, colour=label.call.gt)) + 
		scale_colour_manual(values = colour.call.gt) +
		scale_x_continuous(expand = c(0.005, 0.005)) +
		#coord_cartesian(ylim = c(-7.5, 102.5)) +
		theme_classic() +
		theme(panel.grid.major = element_line(colour="grey90"),
					panel.grid.major.x = element_blank(),
					strip.text = element_text(face = "bold"),
					legend.title = element_blank(),
					strip.background = element_rect(fill = NA, colour = NA),
					panel.border = element_rect(fill = NA, colour = "black"),
					panel.margin = unit(0.25, units = "inches"),
					panel.margin.y = unit(0.1, units = "inches"),
					legend.position = "bottom") +
		xlab("Chromosome position (Mb)") +
		ylab("Relative error proportion (%)")
}



gg = plot.prop.together(d.prop.together)
ggsave(sprintf("_plot.maf.prop_together.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)

gg = plot.prop.together(d.prop.together[d.prop.together$error,], 5, 11)
ggsave(sprintf("_plot.maf.prop_together_error.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)

gg = plot.size.together(d.prop.together)
ggsave(sprintf("_plot.maf.size_together.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 6)


gg = plot.prop.bytruegt(d.prop.bytruegt)
ggsave(sprintf("_plot.maf.prop_bytruegt.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)

gg = plot.prop.bytruegt(d.prop.bytruegt[d.prop.bytruegt$error, ], 8, 17)
ggsave(sprintf("_plot.maf.prop_bytruegt_error.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)

gg = plot.size.bytruegt(d.prop.bytruegt)
ggsave(sprintf("_plot.maf.size_bytruegt.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 6)


gg = plot.prop.bytruegt.lines(d.prop.bytruegt)
ggsave(sprintf("_plot.maf.lines.prop_bytruegt.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)

gg = plot.prop.bytruegt.lines(d.prop.bytruegt[d.prop.bytruegt$error, ], 8, 17)
ggsave(sprintf("_plot.maf.lines.prop_bytruegt_error.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 12)


gg = plot.pos.together(d.pos.together)
ggsave(sprintf("_plot.pos.together.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 18)

gg = plot.pos.bytruegt(d.pos.bytruegt)
ggsave(sprintf("_plot.pos.bytruegt.%s.typed_%s.pdf", named, typed), gg, width = 12, height = 18)


cat("OK\n")



cat("Stats ... ")

gt = as.character(0:2)
m = matrix(0, nrow = 3, ncol = 3, dimnames = list(gt, gt))
for (tgt in gt) {
	for (cgt in gt) {
		m[tgt, cgt] = length(which(p$true.gt == as.numeric(tgt) & p$call.gt == as.numeric(cgt)))
	}
}
m = m / apply(m, 1, sum)

save(m, file = sprintf("_stat.prop_matrix.%s.typed_%s.RData", named, typed))

cat("OK\n")






d00 = data.frame(x=g0, y = g0, true="0", tag="0 as 0")
d01 = data.frame(x=g0, y = g1, true="0", tag="0 as 1")
d02 = data.frame(x=g0, y = g2, true="0", tag="0 as 2")

d10 = data.frame(x=g1, y = g0, true="1", tag="1 as 0")
d11 = data.frame(x=g1, y = g1, true="1", tag="1 as 1")
d12 = data.frame(x=g1, y = g2, true="1", tag="1 as 2")

d20 = data.frame(x=g2, y = g0, true="2", tag="2 as 0")
d21 = data.frame(x=g2, y = g1, true="2", tag="2 as 1")
d22 = data.frame(x=g2, y = g2, true="2", tag="2 as 2")


d = rbind(d00, d01, d02,  d10, d11, d12,  d20, d21, d22)

ggplot(data = d) + facet_wrap(~true) + geom_point(aes(x=x, y=y, colour=tag))




