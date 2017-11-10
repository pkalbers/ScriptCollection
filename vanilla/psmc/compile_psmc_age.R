

library(data.table)


tt = fread("../../vanilla.times.txt", header = T, stringsAsFactors = F)


d = NULL
t = NULL

files = dir(pattern = "^psmc_rvage_.+\\.RData$")

for (file in files) {
	
	load(file)
	
	tag = sub("^.+_pack\\.([0-9a-zA-Z]{16})\\.txt\\..+\\.RData$", "\\1", file)
	nnr = sub("^psmc_rvage_([A-Z]{2})_pack\\..+", "\\1", file)
	
	if (is.null(t)) {
		t = times
	} else {
		if (!identical(t, times)) {
			stop("!!!")
		}
	}
	
	names(psmc) = sprintf("%s %s %s", names(psmc), tag, nnr)
	
	d = c(d, psmc)
	
}


rt = t + ( diff(c(t, 15)) / 2 )


p = data.frame(key = names(d), 
							 mid = as.numeric(sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\1", names(d))),
							 pos = as.numeric(sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\2", names(d))),
							 id0 = as.numeric(sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\3", names(d))),
							 id1 = as.numeric(sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\4", names(d))),
							 shr = as.numeric(sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\5", names(d))),
							 tag = (sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\6", names(d))),
							 nnr = (sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\7", names(d))),
							 stringsAsFactors = F)
p = data.table(p)








exclude.threshold = function(x) {
	a = which(x$shr == 1)
	b = which(x$shr == 0)
	
	ti = exp(seq(log(1e-8), log(40), length.out = 1024))
	
	# q25 = NULL
	# for (t in ti) {
	# 	if (min(x$q25) > t) next
	# 	if (max(x$q25) < t) next
	# 	q25 = rbind(q25, data.table(t = t, a = length(which(x$q25[a] > t)), b = length(which(x$q25[b] < t))))
	# }
	q50 = NULL
	for (t in ti) {
		if (min(x$q50) > t) next
		if (max(x$q50) < t) next
		q50 = rbind(q50, data.table(t = t, a = length(which(x$q50[a] > t)), b = length(which(x$q50[b] < t))))
	}
	# q75 = NULL
	# for (t in ti) {
	# 	if (min(x$q75) > t) next
	# 	if (max(x$q75) < t) next
	# 	q75 = rbind(q75, data.table(t = t, a = length(which(x$q75[a] > t)), b = length(which(x$q75[b] < t))))
	# }
	
	#q25$n = q25$a / length(a) + q25$b / length(b)
	q50$n = q50$a / length(a) + q50$b / length(b)
	#q75$n = q75$a / length(a) + q75$b / length(b)
	
	# list(q25 = q25$t[order(q25$n)[1]],
	# 		 q50 = q50$t[order(q50$n)[1]],
	# 		 q75 = q75$t[order(q75$n)[1]])
	q50$t[order(q50$n)[1]]
}

exclude.pairs = function(stat, y) {
	stat$pass = 1
	
	th = exclude.threshold(stat)
	
	a = which(stat$shr == 1)
	b = which(stat$shr == 0)
	
	if (length(a) > 1) {
		aa = which(stat$q50[a] > th)
		if (length(aa) > 0) {
			xa = stat$key[a[aa]]
			if (length(aa) > length(a) / 2) {
				xa = xa[ order(stat$q50[a[aa]], decreasing = T) ]
				xa = xa[ 1:(floor(length(a) / 2))]
			}
			xa = which(stat$key %in% xa)
			stat$pass[xa] = 0
		}
	}
	
	if (length(b) > 1) {
		bb = which(stat$q50[b] < th)
		if (length(bb) > 0) {
			xb = stat$key[b[bb]]
			if (length(bb) > length(b) / 2) {
				xb = xb[ order(stat$q50[b[bb]]) ]
				xb = xb[ 1:(floor(length(b) / 2))]
			}
			xb = which(stat$key %in% xb)
			stat$pass[xb] = 0
		}
	}
	
	k = which(stat$pass == 1)
	k = stat$key[k]
	
	list(stat=stat, y=y[k])
}




age = list()

m = split(p, sprintf("%d %s %s", p$mid, p$tag, p$nnr))

for (mm in names(m)) {
	
	for (adj in c(0, 1)) {
		
		x = m[[mm]]
		x = x[order(x$key),]
		y = d[x$key]
		
		stat = lapply(y, function(p, rt) {
			cum = cumsum(p)
			data.table(mean = sum(rt * p),
								 mode = rt[ which.max(p) ],
								 q25 = rt[ which.min(abs(cum - 0.25)) ],
								 q50 = rt[ which.min(abs(cum - 0.50)) ],
								 q75 = rt[ which.min(abs(cum - 0.75)) ])
		}, rt)
		stat = rbindlist(stat)
		stat = cbind(x, stat)
		stat$pass = 1
		
		if (adj == 1) {
			sy = exclude.pairs(stat, y)
			stat = sy$stat
			y = sy$y
		}
		
		y = lapply(y, function(p) {
			p = p / sum(p)
			cumsum(p)
		})
		
		y = lapply(y, function(p) {
			if (any(p > 1)) {
				w = which(p > 1)
				p[w] = 1
			}
			p
		})
		
		z = stat$key[which(stat$shr == 0 & stat$pass == 1)]
		
		y[z] = lapply(y[z], function(p) {
			1 - p
		})
		
		y = lapply(y, function(p) {
			log(p)
		})
		
		if (any(is.na(unlist(y)))) stop("!!!!")
		
		cle = 0
		
		for (k in 1:length(y)) {
			cle = cle + y[[k]]
		}
		
		seq = exp(cle)
		seq = seq / sum(seq)
		cum = cumsum(seq)
		
		tmp = list()
		
		tmp$stat = stat
		
		tmp$cle = data.table(t = t, p = exp(cle - max(cle)))
		
		y = lapply(y, function(p, t) {
			data.table(t = t, p = exp(p))
		}, t)
		
		for (k in names(y)) {
			y[[k]]$pair = sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\3 \\4", k)
			s = sub("^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-1]+) ([0-9a-zA-Z]{16}) ([A-Z]{2})$", "\\5", k)
			y[[k]]$cord = ""
			if (s == "1") { y[[k]]$cord = "Concordant" }
			if (s == "0") { y[[k]]$cord = "Discordant" }
		}
		
		tmp$ccf = rbindlist(y)
		
		tmp$tru = tt[ as.numeric(stat$mid[1])+1 , ]
		
		tmp$key = data.table(mid = stat$mid[1], tag = stat$tag[1], nnr = stat$nnr[1], adj = adj)
		
		tmp$est = data.table(mean = sum(rt * seq),
												 mode = rt[ which.max(tmp$cle$p) ],
												 q025 = rt[ which.min(abs(cum - 0.025)) ],
												 q50 = rt[ which.min(abs(cum - 0.50)) ],
												 q975 = rt[ which.min(abs(cum - 0.975)) ])
		
		age[[ sprintf("%s %d", mm, adj) ]] = tmp
	}
}



save(age, file = "_result_psmc.RData")


stop()
#####


library(ggplot2)
library(ggthemes)




a = age[[1]]

nnr = a$stat$nnr[1]
mid = a$stat$mid[1]



gg = ggplot(a$ccf) +
	geom_step(aes(x = (t + 1e-8) * 20000, y=p, color = cord, group = pair), alpha = 1) +
	geom_area(data = data.table(x = c(1e-8, 15) * 20000, y = c(1, 1) + 0.1), aes(x, y), fill = "white", alpha = 0.5) + 
	geom_step(data = a$cle, aes(x = (t + 1e-8) * 20000, y=p), alpha = 3/4, size = 2, color = "white") +
	geom_step(data = a$cle, aes(x = (t + 1e-8) * 20000, y=p), alpha = 3/4, size = 1) +
	geom_vline(xintercept = c(a$tru$node.time, a$tru$parent.time), color = "blue", alpha = 2/3, size = 1) +
	geom_vline(xintercept = c(a$est$q025,a$est$q975) * 20000, linetype = "11", alpha = 2/3, size = 1) +
	scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = format(c(1, 10, 100, 1000, 10000, 100000), scientific = F, big.mark = ',')) +
	scale_colour_manual(values = c("Concordant" = "forestgreen", "Discordant" = "chocolate")) +
	coord_cartesian(xlim = c(0.002, 10) * 20000, ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(panel.border = element_rect(fill = NA, colour = "black", size = 2/3), 
				legend.position = "none") +
	xlab("Time (generations)") + ylab("CCF")
gg

ggsave(gg, filename = sprintf("_plot_psmc_full.%d.%s.pdf", mid, nnr), width = 297, height = 210, units = "mm")



#######################







library(ggplot2)
library(ggthemes)

Ne = 2 * 10000



marker = fread("../../vanilla.marker.txt", header = T, stringsAsFactors = F)




for (ttt in sample(names(age), 25)) {
	
	a = age[[ttt]]
	
	
	tag = a$stat$tag[1]
	nnr = a$stat$nnr[1]
	mid = a$stat$mid[1]
	
	f.ccf = dir(pattern = sprintf(".+_%s_.+\\.%s\\..+\\.ccf\\.txt$", nnr, tag), path = "../rvage", full.names = T)
	f.cle = dir(pattern = sprintf(".+_%s_.+\\.%s\\..+\\.cle\\.txt$", nnr, tag), path = "../rvage", full.names = T)
	
	d.ccf = fread(f.ccf, header = T, stringsAsFactors = F)
	d.cle = fread(f.cle, header = T, stringsAsFactors = F)
	
	d.ccf = d.ccf[ which(d.ccf$MarkerID == mid) ]
	d.cle = d.cle[ which(d.cle$MarkerID == mid) ]
	
	
	pair = rbind(data.table(MarkerID = d.ccf$MarkerID,
													Clock = d.ccf$Clock,
													Pair = sprintf("%d %d", (d.ccf$SampleID0 * 2) + d.ccf$Chr0, (d.ccf$SampleID1 * 2) + d.ccf$Chr1),
													Shared = d.ccf$Shared, 
													q25 = d.ccf$q25,
													q50 = d.ccf$q50,
													q75 = d.ccf$q75), 
							 data.table(MarkerID = a$stat$mid,
							 					 Clock = "PSMC",
							 					 Pair = sprintf("%d %d", a$stat$id0, a$stat$id1),
							 					 Shared = a$stat$shr, 
							 					 q25 = a$stat$q25,
							 					 q50 = a$stat$q50,
							 					 q75 = a$stat$q75))
	
	range(table(pair$Pair))
	
	site = rbind(data.table(MarkerID = d.cle$MarkerID,
													Clock = d.cle$Clock,
													Adjusted = d.cle$Adjusted,
													Fk = d.cle$Fk,
													PostMean = d.cle$PostMean,
													PostMode = d.cle$PostMode,
													PostMedian = d.cle$PostMedian,
													PostCI025 = d.cle$PostCI025,
													PostCI975 = d.cle$PostCI975),
							 data.table(MarkerID = d.cle$MarkerID[1],
							 					 Clock = "PSMC",
							 					 Adjusted = 0,
							 					 Fk = d.cle$Fk[1],
							 					 PostMean = a$est$mean,
							 					 PostMode = a$est$mode,
							 					 PostMedian = a$est$q50,
							 					 PostCI025 = a$est$q025,
							 					 PostCI975 = a$est$q975))
	
	
	tru = data.table(x = c(-1, 3), ymin = a$tru$node.time, ymax = a$tru$parent.time)
	
	
	
	
	x = pair
	y = site
	
	near = "Nearest neighbour"
	if (nnr == "RD") near = "Random"
	
	title = sprintf("Allele count: %d     Pair selection: %s", a$tru$leaves, near)
	
	
	
	y$fltr = "Raw"
	y$fltr[which(y$Adjusted == 1)] = "Adjusted"
	y$fltr = factor(y$fltr, levels = c("Raw", "Adjusted"), labels = c("Raw", "Adjusted"), ordered = T)
	
	x = rbindlist(lapply(split(x, x$Clock), function(x) {
		x = x[order(abs(x$Shared-1), x$q50),]
		
		a = which(x$Shared == 1)
		b = which(x$Shared == 0)
		
		x$i = NA
		x$i[a] = ((1:length(a)) / (length(a)+1)) - 0.06
		x$i[b] = 1 + ((1:length(b)) / (length(b)+1)) + 0.06
		
		x
	}))
	
	
	x$cord = NA
	x$cord[ which(x$Shared == 1) ] = "Concordant pairs"
	x$cord[ which(x$Shared == 0) ] = "Discordant pairs"
	
	x$Clock[ which(x$Clock == "M") ] = "(a) Mutation clock"
	x$Clock[ which(x$Clock == "R") ] = "(b) Recombination clock"
	x$Clock[ which(x$Clock == "C") ] = "(c) Combined clock"
	x$Clock[ which(x$Clock == "PSMC") ] = "(d) PSMC"
	
	y$Clock[ which(y$Clock == "M") ] = "(a) Mutation clock"
	y$Clock[ which(y$Clock == "R") ] = "(b) Recombination clock"
	y$Clock[ which(y$Clock == "C") ] = "(c) Combined clock"
	y$Clock[ which(y$Clock == "PSMC") ] = "(d) PSMC"
	
	
	
	
	
	gg = ggplot(x) +
		facet_grid(.~Clock) +
		geom_rect(data = data.table(xmin = 0.95, xmax = 1.05, ymin = 1e-8, ymax = 1e8), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", alpha = 1/2) +
		geom_ribbon(data = tru, aes(x=x, ymin=ymin, ymax=ymax), fill = "grey70", color = "grey50", alpha = 1/2, size = 1/3, linetype = "12", show.legend = F) +
		geom_ribbon(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, fill = cord), alpha = 2/3, show.legend = F) +
		geom_line(aes(x=i, y=q50 * Ne, color = cord), size = 1/2, alpha = 1, show.legend = F) +
		#geom_linerange(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, color = concord), alpha = 1/2, size = 1/3, show.legend = F) +
		#geom_point(aes(x=i, y=q50 * Ne, color = concord), shape = 16, alpha = 1, size = 1/2, show.legend = F) +
		#geom_vline(data = q, aes(xintercept = x), size = 1/3, alpha = 1/2) +
		#geom_rect(data = y, aes(xmin = 0.955, xmax = 1.045, ymin=PostCI025 * Ne, ymax=PostCI975 * Ne, fill = fltr), alpha = 1/2) +
		#geom_segment(data = y, aes(x = 0.956, xend = 1.044, y = PostMode * Ne, yend = PostMode * Ne, linetype = fltr), alpha = 1/2, size = 1/2) +
		geom_pointrange(data = y, aes(x = 1, y = PostMode * Ne, ymin=PostCI025 * Ne, ymax=PostCI975 * Ne, shape = fltr), alpha = 1, position = position_dodge(width = 0.1), fatten = 3) +
		#geom_segment(data = q, aes(x = x0, xend = x1, y = y0, yend = y1 * Ne), alpha = 1/10, size = 1/4) +
		scale_x_continuous(breaks = c(0.45, 1, 1.55), labels = c("Concordant\npairs", "Age\nestimate", "Discordant\npairs")) +
		scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ','), minor_breaks = c(seq(1,10,by=1), seq(10,100,by=10), seq(100,1000,by=100), seq(1000,10000,by=1000), seq(10000,100000,by=10000))) +
		scale_colour_manual(values = c("Concordant pairs" = "forestgreen", "Discordant pairs" = "chocolate")) +
		scale_fill_manual(values = c("Concordant pairs" = "forestgreen", "Discordant pairs" = "chocolate")) +
		scale_shape_manual(values = c("Raw" = 1, "Adjusted" = 19)) +
		scale_linetype_manual(values = c("Raw" = "11", "Adjusted" = "solid")) +
		coord_cartesian(ylim = c(0.8, 1.2e5), xlim = c(-0.075, 2.075), expand = F) +
		theme_few() +
		theme(aspect.ratio = 16/6,
					panel.background = element_rect(fill = "grey85"),
					panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
					strip.text.x = element_text(face = "bold", hjust = 0, size = 12),
					axis.title.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.text.x = element_text(size = 8),
					axis.title.y = element_text(margin = margin(1,-1,1,1, "mm")),
					plot.title = element_text(size = 10, hjust = 1),
					plot.margin = margin(1, 1, 1, 1, "cm"),
					axis.line = element_blank(),
					legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
					legend.key.size = unit(4.25, "mm"), 
					legend.text = element_text(size = 8),
					legend.justification = c(1, 0), legend.position = c(0.995, 1-0.991), 
					legend.title = element_blank(),
					legend.margin = margin(-0.25,2,1,1, "mm"),
					panel.grid.major.y = element_line(colour = "grey99", size = 1/2),
					panel.grid.minor.y = element_line(colour = "grey99", size = 1/3),
					panel.grid.major.x = element_blank(),
					panel.grid.minor.x = element_blank()) +
		ylab("Time (generations)") +
		ggtitle(title)
	gg
	
	ggsave(gg, filename = sprintf("_plot_psmc_rvage.%d.%s.pdf", x$MarkerID[1], nnr), width = 297, height = 210, units = "mm")
	
	
	
}






########################



library(data.table)
library(ggplot2)
library(ggthemes)



Ne = 2 * 10000

marker = fread("../../vanilla.marker.txt", header = T, stringsAsFactors = F)
times = fread("../../vanilla.times.txt", header = T, stringsAsFactors = F)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a)
	b = log10(b)
	sqrt(sum((a - b)^2) / n)
}

mle = function(a, b) {
	a = log10(a)
	b = log10(b)
	10^(mean(abs(a - b)))
}

se <- function(x) sqrt(var(x)/length(x))

expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}






q = lapply(age, function(a) {
	cbind(a$key, a$est)
})
q = rbindlist(q)


files = dir(pattern = "^.+\\.cle\\.txt$", path = "../rvage", full.names = T)

d = list()

for (file in files) {
	tmp = fread(file, header = T, stringsAsFactors = F)
	tmp$tag = sub("^.+_pack\\.([0-9a-zA-Z]{16})\\.txt\\..+$", "\\1", basename(file))
	tmp$nnr = sub("^rvage_([A-Z]{2})_pack\\..+", "\\1", basename(file))
	d[[file]] = tmp
}

d = rbindlist(d)


del = which(d$Adjusted == 1)
d = d[-del,]

del = which(q$adj == 1)
q = q[-del,]


p = data.table(MarkerID = q$mid,
							 Clock = "PSMC",
							 Adjusted = q$adj,
							 Fk = NA,
							 N_Shared = NA, 
							 N_Others = NA, 
							 PostMean = q$mean,
							 PostMode = q$mode,
							 PostMedian = q$q50,
							 PostCI025 = q$q025,
							 PostCI975 = q$q975,
							 Robust = NA,
							 Lower = NA,
							 Upper = NA,
							 tag = q$tag,
							 nnr = q$nnr)

p = rbind(d, p)


p = cbind(p, times[p$MarkerID + 1, ])


p$Clock[ which(p$Clock == "M") ] = "Mutation clock"
p$Clock[ which(p$Clock == "R") ] = "Recombination clock"
p$Clock[ which(p$Clock == "C") ] = "Combined clock"
p$Clock[ which(p$Clock == "PSMC") ] = "PSMC"

p$nnr[ which(p$nnr == "NN") ] = "Nearest neighbours"
p$nnr[ which(p$nnr == "RD") ] = "Random pairs"


p$Clock = factor(p$Clock, levels = c("Mutation clock", "Recombination clock", "Combined clock", "PSMC"), labels = c("Mutation clock", "Recombination clock", "Combined clock", "PSMC"), ordered = T)


del = which(p$nnr == "Nearest neighbours")
p = p[-del,]


p$mid = exp(log(p$node.time + 1e-8) + ((log(p$parent.time + 1e-8) - log(p$node.time + 1e-8)) / 2))
p$x = p$mid # p$node.time

mx = c(0.5, 2e5)

tm = rbind(data.table(x = 10^(0:6), y = 1e-6, xend = 10^(0:6), yend = 0.7),
					 data.table(y = 10^(0:6), x = 1e-6, yend = 10^(0:6), xend = 0.7))

tn = c((2:9), 10*(2:9), 100*(2:9), 1000*(2:9), 10000*(2:9))
tn = rbind(data.table(x = tn, xend = tn, y = 1e-6, yend = 0.615), 
					 data.table(y = tn, yend = tn, x = 1e-6, xend = 0.615))

tx = rbind(data.table(x = 1e-6, xend = 1e6, y = 1e-6, yend = 0.7),
					 data.table(x = 1e-6, xend = 0.7, y = 1e-6, yend = 1e6))


w = split(p, sprintf("%s %s", p$Clock, p$nnr))
w = lapply(w, function(x) {
	x$rankc = cor(x$x, x$PostMode, method = "spearman")
	x$rmsle = rmsle(x$x + 1e-32, x$PostMode * Ne + 1e-32)
	x$props = length(which(x$PostMode * Ne >= x$node.time & x$PostMode * Ne <= x$parent.time)) / nrow(x)
	x[1,]
})
w = rbindlist(w)

w$rankc.txt = sprintf("paste('Spearman\\'s ', r[S],' = %.3f ')", w$rankc)
w$rmsle.txt = sprintf("RMSLE = %.3f ", w$rmsle)
w$props.txt = sprintf("Inside interval = %.1f%%", w$props*100)


scat = ggplot(p) +
	#facet_grid(.~Clock) +
	facet_wrap(~Clock, nrow = 2, scales = "free") +
	geom_raster(aes(x, PostMode * Ne, fill = (..density.. / max(..density..)) * 100 ), position = "identity", stat = "bin2d", binwidth = c(0.075, 0.075)) +
	geom_tile(aes(x, PostMode * Ne, fill = (..density.. / max(..density..)) * 100), position = "identity", stat = "bin2d", binwidth = c(0.075, 0.075)) +
	geom_smooth(aes(x, node.time + 1e-8), method = "lm", size = 1/2, color = "grey30", fullrange = F, se = T) +
	#geom_smooth(aes(x, mid), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	geom_smooth(aes(x, parent.time + 1e-8), method = "lm", size = 1/2, color = "grey30", fullrange = F, se = T) +
	geom_text(data = w, aes(x = 150000, y = 4.8, label = rankc.txt), size = 3, hjust = 1, parse = T) +
	geom_text(data = w, aes(x = 150000, y = 2.4, label = rmsle.txt), size = 3, hjust = 1) +
	geom_text(data = w, aes(x = 150000, y = 1.2, label = props.txt), size = 3, hjust = 1) +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = F) +
	geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", se = T, fullrange = F) +
	geom_rect(data = tx, aes(xmin=x, xmax=xend, ymin=y, ymax=yend), fill="grey80") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), size = 1/2, color="black") +
	geom_segment(data = tn, aes(x=x, xend=xend, y=y, yend=yend), size = 1/4, color="black") +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), na.value = "grey90", trans = "log10", limits = c(1, 100), breaks = c(1, 10, 100), labels = c("<1%", "10%", "100%")) +
	scale_color_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), na.value = "grey90", trans = "log10", limits = c(1, 100), breaks = c(1, 10, 100), labels = c("<1%", "10%", "100%")) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 9),
				#axis.ticks = element_blank(),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", size = 11),
				strip.text.y = element_text(size = 10),
				legend.justification=c(0,1), legend.position=c(1-0.98,0.99),
				legend.background = element_rect(fill = NA, colour = NA, size = NA),
				#legend.position = "bottom",
				plot.title = element_text(face = "bold", size = 12),
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 7),
				legend.title = element_blank()) +
	ylab("Inferred age (generations)") + xlab("True age (generations)") #+
	#ggtitle(p$data[1])
#scat$theme$panel.spacing.x = unit(c(0.25, 0.25, 0.75), "cm")
scat

ggsave(scat, filename = ("___plot_psmc_rvage.pdf"), height = 8, width = 8)

ggsave(scat, filename = ("__plot_psmc_rvage.1000.cle.adj.pdf"), height = 7.5, width = 13)




pp = p

pp$grp = "" #cut(pp$parent.time, breaks = c(0, 1000, 1e6), include.lowest = T)

w0 = which(pp$parent.time <= 1000)
w1 = which(pp$node.time >= 1000)
wi = which(pp$node.time < 1000 & pp$parent.time > 1000)

pp$grp[w0] = "a"
pp$grp[wi] = "b"
pp$grp[w1] = "c"


x = by(pp, list(pp$Clock, pp$grp), nrow)
(array(x, dim = dim(x), dimnames = dimnames(x)))

x = by(pp, list(pp$Clock, pp$grp), function(x) cor(x$x, x$PostMode, method = "spearman"))
(array(x, dim = dim(x), dimnames = dimnames(x)))

x = by(pp, list(pp$Clock, pp$grp), function(x) rmsle(x$x, x$PostMode * Ne))
(array(x, dim = dim(x), dimnames = dimnames(x)))

x = by(pp, list(pp$Clock, pp$grp), function(x) length(which(x$PostMode * Ne >= x$node.time & x$PostMode * Ne <= x$parent.time)) / nrow(x) * 100)
(array(x, dim = dim(x), dimnames = dimnames(x)))










######################## 
# True TRMCA
######################## 


library(data.table)
library(ggplot2)
library(ggthemes)



Ne = 2 * 10000

marker = fread("../../vanilla.marker.txt", header = T, stringsAsFactors = F)
times = fread("../../vanilla.times.txt", header = T, stringsAsFactors = F)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a)
	b = log10(b)
	sqrt(sum((a - b)^2) / n)
}




tru = fread("true_ibd_pairs.txt", header = T, stringsAsFactors = F)

tru$id0 = (tru$SampleID0 * 2) + tru$Chr0
tru$id1 = (tru$SampleID1 * 2) + tru$Chr1

tru$key = sprintf("%d %d %d %d %d", tru$MarkerID, tru$SampleID0, tru$Chr0, tru$SampleID1, tru$Chr1)
tru$tag = sprintf("%d %d %d", tru$MarkerID, tru$id0, tru$id1)



load("../psmc/_result_psmc.RData")

q = lapply(age, function(a) {
	a$stat
})
q = rbindlist(q)

del = which(q$nnr == "NN")
q = q[-del,]

q$tag = sprintf("%d %d %d", q$mid, q$id0, q$id1)

q = q[sample(1:nrow(q)),]

k = intersect(q$tag, tru$tag)
q = q[match(k, q$tag),]
t = tru[match(k, tru$tag),]
q$tmrca = t$TMRCA




files = dir(pattern = "^.+\\.ccf\\.txt$", path = "../rvage", full.names = T)

d = list()
for (file in files) {
	tmp = fread(file, header = T, stringsAsFactors = F)
	tmp$tag = sub("^.+_pack\\.([0-9a-zA-Z]{16})\\.txt\\..+$", "\\1", basename(file))
	tmp$nnr = sub("^rvage_([A-Z]{2})_pack\\..+", "\\1", basename(file))
	d[[file]] = tmp
}
d = rbindlist(d)

del = which(d$nnr == "NN")
d = d[-del,]

d$key = sprintf("%d %d %d %d %d", d$MarkerID, d$SampleID0, d$Chr0, d$SampleID1, d$Chr1)
d$tag = sprintf("%d %d %d", d$MarkerID, d$SampleID0 * 2 + d$Chr0, d$SampleID1 * 2 + d$Chr1)


d = split(d, d$Clock)
d = lapply(d, function(d) {
	d = d[sample(1:nrow(d)),]
	w = which(duplicated(d$tag))
	if (length(w) > 0) {
		d = d[-w,]
	}
	d
})
d = rbindlist(d)


k = intersect(Reduce(intersect, lapply(split(d, d$Clock), function(x) x$tag)), q$tag)
q = q[which(q$tag %in% k), ]
d = d[which(d$tag %in% k), ]


d$mode = (d$Shape - 1) / d$Rate
d$mean = (d$Shape) / d$Rate
d$median = d$q50
x = which(d$mode == 0)
if (length(x) > 0) {
	d$mode[x] = d$mode[x] + (0.16/Ne)
}


d = split(d, d$Clock)
d = lapply(d, function(d, t) {
	k = intersect(d$key, t$key)
	d = d[match(k, d$key),]
	t = t[match(k, t$key),]
	d$tmrca = t$TMRCA
	d
}, tru)
d = rbindlist(d)



# del = which(d$Pass == 0)
# d = d[-del,]
# 
# del = which(q$pass == 0)
# q = q[-del,]



p = rbind(data.table(mid = d$MarkerID, clk = d$Clock, est = d$median, tru = d$tmrca, shr = d$Shared, tag = d$tag), 
					data.table(mid = q$mid, clk = "PSMC", est = q$q50, tru = q$tmrca, shr = q$shr, tag = q$tag))

range(table(p$tag))
nrow(p)/4

p$clk[ which(p$clk == "M") ] = "(a) Mutation clock"
p$clk[ which(p$clk == "R") ] = "(b) Recombination clock"
p$clk[ which(p$clk == "C") ] = "(c) Combined clock"
p$clk[ which(p$clk == "PSMC") ] = "(d) PSMC"

p$shr[ which(p$shr == 1) ] = "Concordant pairs"
p$shr[ which(p$shr == 0) ] = "Discordant pairs"


p$clk = factor(p$clk, levels = sort(unique(p$clk)), labels = sub("^.{4}(.+)$", "\\1", sort(unique(p$clk))), ordered = T)


Ne = 20000




mx = c(0.1, 2e5)

tm = rbind(data.table(x = 10^(0:6), y = 1e-6, xend = 10^(0:6), yend = 0.15),
					 data.table(y = 10^(0:6), x = 1e-6, yend = 10^(0:6), xend = 0.15))

tn = c(0.1*(2:9), (2:9), 10*(2:9), 100*(2:9), 1000*(2:9), 10000*(2:9))
tn = rbind(data.table(x = tn, xend = tn, y = 1e-6, yend = sqrt((mx[1]*1.1) * (0.15))), 
					 data.table(y = tn, yend = tn, x = 1e-6, xend = sqrt((mx[1]*1.1) * (0.15))))

tx = rbind(data.table(x = 1e-6, xend = 1e6, y = 1e-6, yend = 0.15),
					 data.table(x = 1e-6, yend = 1e6, y = 1e-6, xend = 0.15))


x = which(p$tru < 0.2)
if (length(x) > 0) {
	p$tru[x] = 0.16
}


pp = p
del = which(pp$est*Ne < 0.2)
if (length(del) > 0) {
	pp = pp[-del,]
}
del = which(pp$tru < 0.2)
if (length(del) > 0) {
	pp = pp[-del,]
}

w = split(pp, sprintf("%s %s", pp$clk, pp$shr))
w = lapply(w, function(x) {
	x$rankc = cor(x$tru, x$est, method = "spearman")
	x$rmsle = rmsle(x$tru, x$est * Ne)
	x[1,]
})
w = rbindlist(w)

w$rankc.txt = sprintf("paste('Spearman\\'s ', r[S],' = %.3f ')", w$rankc)
w$rmsle.txt = sprintf("RMSLE = %.3f ", w$rmsle)


#cd = data.table(shr = sort(unique(p$shr)), lab = c("Concordant", "Discordant"), x = sqrt(prod(mx)), y = 1e5 * 0.9)


scat = ggplot(pp) +
	facet_grid(clk~shr) +
	geom_raster(aes(tru, est * Ne, fill = (..density.. / max(..density..)) ), position = "identity", stat = "bin2d", binwidth = c(0.1, 0.1)) +
	geom_tile(aes(tru, est * Ne, fill = (..density.. / max(..density..)) ), position = "identity", stat = "bin2d", binwidth = c(0.1, 0.1)) +
	geom_text(data = w, aes(x = 150000, y = 1.0, label = rankc.txt), size = 3, hjust = 1, parse = T) +
	geom_text(data = w, aes(x = 150000, y = sqrt(0.15), label = rmsle.txt), size = 3, hjust = 1) +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_smooth(data = pp, aes(tru, est * Ne), method = "glm", color = "black", size = 2/3, se = F, fullrange = F) +
	geom_smooth(data = pp, aes(tru, est * Ne), method = "glm", color = "white", linetype = "22", size = 2/3, fill = "black", se = T, fullrange = F) +
	geom_rect(data = tx, aes(xmin=x, xmax=xend, ymin=y, ymax=yend), fill="grey80") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), size = 1/2, color="black") +
	geom_segment(data = tn, aes(x=x, xend=xend, y=y, yend=yend), size = 1/4, color="black") +
	#geom_label(data = cd, aes(x=x, y=y, label = lab, color = shr), show.legend = F, size = 2.5, label.r = unit(0.25, "lines"), label.padding = unit(0.25, "lines"), alpha = 2/3) + 
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	# scale_fill_gradientn(colours = c("grey90", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), 
	# 										 na.value = "grey90", trans = "log10", limits = c(0.01, 1), breaks = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1")) +
	scale_fill_gradientn(colours = c("grey90", "khaki", "seagreen3", "seagreen", "navy"), 
											 na.value = "grey90", trans = "log10", limits = c(0.01, 1), breaks = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1")) +
	#scale_color_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), na.value = "grey90", trans = "log10", limits = c(1, 100), breaks = c(1, 10, 100), labels = c("<1%", "10%", "100%")) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	scale_colour_manual(values = c("forestgreen", "chocolate")) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 8.5),
				axis.ticks = element_line(color = "black"),
				axis.title.x = element_text(margin = margin(1,0.1,0.1,0.1, unit = "lines")),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				strip.text.x = element_text(size = 12),
				strip.text.y = element_text(face = "bold", size = 12),
				#strip.text.y = element_text(size = 12, margin = margin(0.1,0.1,0.1,1, unit = "lines")),
				legend.justification=c(0,1), legend.position=c(1-0.975,0.995),
				legend.background = element_rect(fill = NA, colour = NA, size = NA),
				#legend.position = "bottom",
				plot.title = element_text(face = "bold", size = 13),
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 7),
				legend.title = element_blank()) +
	ylab("Inferred TMRCA (generations)") + xlab("True TMRCA (generations)")
scat$theme$panel.spacing.y = unit(c(0.25, 0.25, 0.75), "cm")
scat


ggsave(scat, filename = ("__plot_psmc_rvage.1000.ccf.pdf"), height = 12.5, width = 7)

ggsave(scat, filename = ("__plot_psmc_rvage.1000.ccf.pdf"), height = 7, width = 13)
ggsave(scat, filename = ("__plot_psmc_rvage.1000.ccf.png"), height = 7, width = 13)





table(pp$clk)
table(pp$clk, pp$shr)


pp$grp = cut(pp$tru, breaks = c(0, 100, 1000, 10000, 1e5, 1e6), include.lowest = T)


x = by(pp, list(pp$clk, pp$shr), function(x) cor(x$tru, x$est, method = "spearman"))
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(pp, list(pp$clk, pp$shr), function(x) rmsle(x$tru, x$est * Ne))
array(x, dim = dim(x), dimnames = dimnames(x))



x = by(pp, list(pp$grp, pp$clk, pp$shr), nrow)
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(pp, list(pp$grp, pp$clk, pp$shr), function(x) cor(x$tru, x$est, method = "spearman"))
array(x, dim = dim(x), dimnames = dimnames(x))

x = by(pp, list(pp$grp, pp$clk, pp$shr), function(x) rmsle(x$tru, x$est * Ne))
array(x, dim = dim(x), dimnames = dimnames(x))







