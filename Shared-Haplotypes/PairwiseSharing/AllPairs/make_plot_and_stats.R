

library(ggplot2)
library(grid)
library(data.table)


tab = fread("sub2.chr20_f25.pairwise", header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "integer", "integer", "double", "character",
																																															"character", "character", "character", "character", "character", "character"))

tab = as.data.frame(tab)

last = tab$target[nrow(tab)]
tab = tab[-(which(tab$target == last)), ]

par = list()

for (i in grep("\\|", tab[1, ])) {
	if (grepl("\\:", tab[1, i])) next
	tag = names(tab)[i]
	if (tag %in% c("lpos", "rpos")) {
	cat(tag, "\n")
	spt = strsplit(tab[, i], "|", TRUE)
	spt = lapply(spt, as.numeric)
	par[[tag]] = spt
	}
}


# target sharing plot

brk.max = 5

target = sample(unique(tab$target), 1)

key = which(tab$target == target)

other = tab$other[key]

f = tab$f[key]

foc = tab$mpos[key]

beg = sapply(par$lpos[key], function(x, max) rev(sort(x))[max], brk.max)
end = sapply(par$rpos[key], function(x, max) sort(x)[max], brk.max)

brk.l = lapply(par$lpos[key], function(x, max) rev(sort(x))[1:max], brk.max)
brk.r = lapply(par$rpos[key], function(x, max) sort(x)[1:max], brk.max)



d.esh = data.table(rvs = as.character(foc), other = other, f = f, foc = foc, beg = beg, end = end)
d.esh = d.esh[order(d.esh$foc), ]


d.brk = list()
for (i in 1:length(key)) {
	d.brk[[i]] = data.table(rvs = as.character(foc[i]), other = other[i], brk = unique(c(brk.l[[i]], brk.r[[i]])))
}
d.brk = rbindlist(d.brk)


lim = as.character(foc[which(foc >= 6.5e06 & foc <= 7.5e06)])
d.esh = d.esh[which(d.esh$rvs %in% lim), ]
d.brk = d.brk[which(d.brk$rvs %in% lim), ]
d.esh$rvs = as.numeric(d.esh$rvs)
d.brk$rvs = as.numeric(d.brk$rvs)


gg = ggplot(data=d.esh) +
	facet_grid(rvs~., scales = "free_y", space = "free_y") +
	geom_point(aes(x=foc, y=other), colour="white", size=4, alpha=0.5) +
	geom_point(aes(x=foc, y=other), colour="black", size=3) +
	geom_segment(aes(x=beg, xend=end, y=other, yend=other, group=other), size=1, colour="black") + 
	geom_point(data=d.brk, aes(x=brk, y=other), shape=124, size=4, colour="black") +
	geom_vline(xintercept=c(min(d.esh$beg), max(d.esh$end)), size=1.5) +
	theme_classic() +
	theme(axis.ticks.length=unit(0.5, "lines"),
				axis.title.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.line.y=element_blank(),
				strip.text.y=element_blank(),
				strip.background=element_blank(),
				panel.margin.x=unit(0.5, "lines"),
				panel.background=element_rect(fill="grey95"),
				legend.position="none") +
	scale_x_continuous(expand = c(0, 0), 
										 breaks=round(seq(0, 60e06, length.out = 121) / 1e06, 2) * 1e06, 
										 labels=round(seq(0, 60e06, length.out = 121) / 1e06, 2)) +
	#coord_cartesian(ylim=c(6.5e06, 7.5e06)) +
	xlab("Chromosome position (Mb)")

ggsave(gg, filename = sprintf("plot.ESH.%s.png", target), width = 10, height = length(unique(paste(d.esh$rvs, d.esh$other))) / 5, limitsize=FALSE)



# cumsum of f's

ftab = table(tab$f)
ftab = cumsum(ftab)
d = data.table(f = as.numeric(names(ftab)), n = ftab)
d$f = factor(d$f)

ggplot(data=d) +
	geom_bar(aes(x=f, y=n/1e06), stat="identity", fill="grey30") +
	theme_classic() +
	ylab("# variants (millions)") +
	xlab("f") +
	ggtitle("Cumulative number of variants") +
	coord_cartesian(ylim=c(0, 20))



# coverage

# indv = unique(c(tab$target, tab$other))
# 
# pair = combn(indv, 2)
# 
# tab.key = paste(tab$target, tab$other)
# pr0.key = paste(pair[1, ], pair[2, ])
# pr1.key = paste(pair[2, ], pair[1, ])
# 
# i = which(tab.key %in% pr0.key)
# j = which(tab.key %in% pr1.key)
# 
# key = unique(c(i, j))
# 
# 
# tab = tab[key, ]
# par = lapply(par, function(p, k) p[k, ], key)



abs.beg = min(sapply(par$lpos, min))
abs.end = max(sapply(par$rpos, max))
dst = abs.end - abs.beg + 1

beg = sapply(par$lpos, max)
end = sapply(par$rpos, min)

targets = unique(tab$target)
targets = sample(targets, 100)

select = which(tab$target %in% targets)

d = cbind(tab[select, ], beg = beg[select], end = end[select])


splt = split(d, d$target)

cumcov = matrix(0, nrow = length(targets), ncol = length(2:25), dimnames = list(targets, as.character(2:25)))

for (target in names(splt)) {
	
	target.tab = splt[[target]]
	
	x = rep(FALSE, dst)
	
	for (f in 2:25) {
		fi = which(target.tab$f == f)
		
		f.tab = target.tab[fi, ]
		
		for (i in 1:nrow(f.tab)) {
			rng = (f.tab$beg[i] - abs.beg + 1):(f.tab$end[i] - abs.beg + 1)
			x[rng] = TRUE #x[rng] + 1
		}
		
		res = length(which(x))
		
		cat(f, target, sprintf("%.3f", res / dst), "\n")
		
		cumcov[target, as.character(f)] = res
	}
	
	cat("\n")
}

save(cumcov, dst, file = "cumcov.RData")


cc = cumcov / dst

d = NULL
for (f in colnames(cc)) {
	d = rbind(d, data.table(f = as.numeric(f), cov = cc[, f]))
}

ex = 1000
gg = ggplot(data=d, aes(x=f, y=ex^cov)) + 
	geom_hline(yintercept=ex) +
	geom_hline(yintercept=ex^c(0.6, 0.7, 0.8, 0.9), linetype="dashed", colour="grey60") +
	geom_hline(yintercept=ex^c(0.95, 0.96, 0.97, 0.98, 0.99), linetype="dotted", colour="grey80") +
	stat_summary(fun.data=mean_cl_normal, geom='smooth', colour="black",  fill="royalblue3", alpha=0.75) +
	stat_summary(fun.y=mean, na.rm = TRUE, geom='point', shape=21, colour="white", fill="black", size=1.5) +
	#geom_boxplot(aes(x=factor(f), y=len), fill="grey", outlier.size=0) +
	theme_classic() +
	scale_x_continuous(breaks=2:25, labels=as.character(2:25)) + 
	#scale_x_discrete(expand=c(0,0)) +
	scale_y_continuous(breaks=ex^c(0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1), labels = c(0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1) * 100) +
	coord_cartesian(ylim=ex^c(0.3, 1)) +
	xlab("f") +
	ylab("Mean coverage ± 95% CI (%)") +
	ggtitle("Cumulative coverage")
ggsave(gg, filename = "plot.1KG_cumulative_coverage.pdf", width = 9+1, height = 16, scale = 1/2.5)



# median length of shared segments

beg = sapply(par$lpos, max)
end = sapply(par$rpos, min)

d = data.table(tab, len = end - beg + 1)
#d = d[-(which(d$len < 0))]


gg = ggplot(data=d, aes(x=f, y=len / 1e06)) + 
	#stat_summary(fun.y=median, na.rm = TRUE, geom='point') +
	#stat_summary(fun.data=mean_cl_normal, geom='smooth', color='black') +
	geom_boxplot(aes(x=factor(f), y=len), fill="grey", outlier.size=0) +
	theme_classic() +
	scale_x_continuous(breaks=2:25, labels=as.character(2:25)) + 
	#scale_x_discrete(expand=c(0,0)) +
	#scale_y_log10(expand=c(0, 0), limits=c(0, Inf)) +
	coord_cartesian(ylim=c(0, 1)) +
	xlab("f") +
	ylab("Median physical distance (Mb)") +
	ggtitle("Lengths of detected shared haplotype segments")
ggsave(gg, filename = "plot.1KG_segment_length.pdf", width = 9+1, height = 16, scale = 1/2.5)




# rate of RVs per 10Mb

beg = min(sapply(par$lpos, min))
end = max(sapply(par$rpos, max))

targets = unique(tab$target)
targets = targets[1:100]

x = list()

for (t in 1:100) {
	target = targets[t]
	
	key = which(tab$target == target)
	
	f = tab$f[key]
	p = tab$mpos[key]
	
	z = which(duplicated(p))
	f = f[-z]
	p = p[-z]
	
	b = cut(p, breaks = seq(0, 100e6, by=1e6), include.lowest = TRUE)
	
	splt = split(f, b)
	
	x[[t]] = lapply(splt, table)
}


R = list()

for (t in 1:100) {
	r = 2:25
	flv = as.character(r)
	names(r) = flv
	
	l = length(x[[t]])
	
	for (i in 1:l) {
		b = x[[t]][[i]]
		
		for (f in flv) {
			if (!is.na(b[f])) {
				r[f] = r[f] + b[f]
			}
		}
	}
	
	r = r /l
	
	R[[t]] = r
}

R = Reduce(rbind, R)


d = list()
for (col in colnames(R)) {
	d[[col]] = data.table(f=as.numeric(col), r=R[, col])
}
d = rbindlist(d)

ggplot(data=d, aes(x=f, y=r)) +
	stat_summary(fun.y=mean, na.rm = TRUE, geom='point') + #, fill=NA, size=4, shape=21) + 
	stat_summary(fun.data=mean_cl_normal, geom='smooth', se=T, color='black') +
	theme_classic() +
	scale_x_continuous(breaks=2:25, labels=as.character(2:25))


CR = t(apply(R, 1, cumsum))

d = list()
for (col in colnames(CR)) {
	d[[col]] = data.table(f=as.numeric(col), r=CR[, col])
}
d = rbindlist(d)

ggplot(data=d, aes(x=f, y=r)) +
	stat_summary(fun.y=mean, na.rm = TRUE, geom='point') + #, fill=NA, size=4, shape=21) + 
	stat_summary(fun.data=mean_cl_normal, geom='smooth', se=T, color='black') +
	theme_classic() +
	scale_x_continuous(breaks=2:25, labels=as.character(2:25)) 




plot( 2:25, colMeans(R) )

plot( 2:25, colMeans( t( apply(R, 1, cumsum) ) ) , xlim=c(2, 25), ylim=c(2, 25))
lines(2:25, 2:25)








