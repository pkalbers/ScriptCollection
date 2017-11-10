


library(data.table)
library(ggplot2)
library(ggthemes)

Ne = 2 * 7300 # 10000

marker = fread("./truH.marker.txt", header = T)
times = fread("~/Research/africa/dev/OutOfAfricaHapMap20.times.txt", header = T)
mut = as.data.table(read.table("~/Research/DetectIBD/OutOfAfricaHapMap20.mutations.txt", header = T))
rec = as.data.table(read.table("~/Research/DetectIBD/OutOfAfricaHapMap20.records.txt", header = T))

rec$c0 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\1", rec$children))
rec$c1 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\2", rec$children))

rec$children = NULL
rec$population = NULL


get.tree = function(tree, rec, pos) {
	out = list()
	for (node in tree$node) {
		sub = which(rec$c0 == node | rec$c1 == node)
		if (length(sub) > 0) {
			int = which(rec$left[sub] <= pos & rec$right[sub] >= pos)
			if (length(int) > 0) {
				out[[as.character(node)]] = get.tree(rec[sub[int],], rec, pos)
			}
		}
	}
	
	if (length(out) == 0) return(tree$node)
	return(unlist(out, use.names = T))
}

get.tmrca = function(pos, h0, h1, rec) {
	loc = which(rec$left <= pos & rec$right >= pos)
	loc = rec[loc, ]
	
	l0 = which(loc$c0 == h0 | loc$c1 == h0)
	l1 = which(loc$c0 == h1 | loc$c1 == h1)
	
	if (length(l0) == 0) return(NA)
	if (length(l1) == 0) return(NA)
	
	z = intersect(l0, l1)
	
	if (length(z) == 0) {
		l0 = loc[l0,]
		l1 = loc[l1,]
		
		z0 = get.tree(l0, loc, pos)
		z1 = get.tree(l1, loc, pos)
		
		if (length(z0) > 1) stop("!!!")
		if (length(z1) > 1) stop("!!!")
		
		z0 = c(as.numeric(strsplit(names(z0), '.', T)[[1]]), z0)
		z1 = c(as.numeric(strsplit(names(z1), '.', T)[[1]]), z1)
		z = intersect(z0, z1)
		z = which(loc$node %in% z)
	}
	
	z = z[ which(loc$left[z] <= pos & loc$right[z] >= pos) ]
	
	min(loc$time[z])
}


exclude = function(x) {
	
	a = which(x$Shared == 1)
	b = which(x$Shared == 0)
	
	tt = exp(seq(log(1e-8), log(40), length.out = 1024))
	
	q25 = NULL
	for (t in tt) {
		if (min(x$q25) > t) next
		if (max(x$q25) < t) next
		q25 = rbind(q25, data.table(t = t, a = length(which(x$q25[a] > t)), b = length(which(x$q25[b] < t))))
	}
	q50 = NULL
	for (t in tt) {
		if (min(x$q50) > t) next
		if (max(x$q50) < t) next
		q50 = rbind(q50, data.table(t = t, a = length(which(x$q50[a] > t)), b = length(which(x$q50[b] < t))))
	}
	q75 = NULL
	for (t in tt) {
		if (min(x$q75) > t) next
		if (max(x$q75) < t) next
		q75 = rbind(q75, data.table(t = t, a = length(which(x$q75[a] > t)), b = length(which(x$q75[b] < t))))
	}
	
	q25$n = q25$a / length(a) + q25$b / length(b)
	q50$n = q50$a / length(a) + q50$b / length(b)
	q75$n = q75$a / length(a) + q75$b / length(b)
	
	list(q25 = q25$t[order(q25$n)[1]],
			 q50 = q50$t[order(q50$n)[1]],
			 q75 = q75$t[order(q75$n)[1]])
}





data = "tru"
expos = "ff"
near = "RD"


none = fread(sprintf("dev_%s_H_%s_100_%s_none.age.sites.txt", data, expos, near), header = T)
prop = fread(sprintf("dev_%s_H_%s_100_%s_prop.age.sites.txt", data, expos, near), header = T)
auto = fread(sprintf("dev_%s_H_%s_100_%s_auto.age.sites.txt", data, expos, near), header = T)

d = fread(sprintf("dev_%s_H_%s_100_%s_none.age.pairs.txt", data, expos, near), header = T)
d = cbind(d, times[d$MarkerID + 1, ])

d$h0 = (d$SampleID0 * 2) + d$Chr0
d$h1 = (d$SampleID1 * 2) + d$Chr1



d = split(d, list(d$MarkerID))


x=d$`192308`

#x = d[[ sample(1:length(d), 1) ]]



x = rbindlist(lapply(split(x, x$Clock), function(x) {
	x$age.none = none$PostMode[ which(none$MarkerID == x$MarkerID[1] & none$Clock == x$Clock[1]) ]
	x$age.prop = prop$PostMode[ which(prop$MarkerID == x$MarkerID[1] & prop$Clock == x$Clock[1]) ]
	x$age.auto = auto$PostMode[ which(auto$MarkerID == x$MarkerID[1] & auto$Clock == x$Clock[1]) ]
	x = x[order(abs(x$Shared-1), x$q50),]
	
	a = which(x$Shared == 1)
	b = which(x$Shared == 0)
	
	x$i = NA
	x$i[a] = (1:length(a)) / (length(a)+1)
	x$i[b] = 1 + ((1:length(b)) / (length(b)+1))
	
	x
}))

x$concord = NA
x$concord[ which(x$Shared == 1) ] = "Concordant"
x$concord[ which(x$Shared == 0) ] = "Discordant"

x$Clock[ which(x$Clock == "M") ] = "(a) Mutation clock"
x$Clock[ which(x$Clock == "R") ] = "(b) Recombination clock"
x$Clock[ which(x$Clock == "C") ] = "(c) Combined clock"

ggplot(x) +
	facet_grid(.~Clock) +
	geom_vline(xintercept = 1, size = 1.5, color = "white") +
	geom_pointrange(aes(x=i, y=q50 * Ne, ymin=q25 * Ne, ymax=q75 * Ne, color = concord), alpha = 0.5, size = 0.5, fatten = 0.5) +
	#geom_point(aes(x=i, y = tmrca), shape = 21) +
	#geom_point(aes(x=i, y = ((Shape)/Rate)*Ne)) +
	geom_hline(yintercept = c(x$node.time[1], x$parent.time[1]), alpha = 0.5) +
	#geom_hline(yintercept = exc$q50 * Ne, linetype = "dashed") +
	geom_hline(aes(yintercept = age.none * Ne), color = "#C0C0C0", size = 1, alpha = 3/4) +
	geom_hline(aes(yintercept = age.prop * Ne), color = "#CA93C9", size = 1, alpha = 3/4) +
	geom_hline(aes(yintercept = age.auto * Ne), color = "#A5A9FF", size = 1, alpha = 3/4) +
	geom_hline(aes(yintercept = age.none * Ne), linetype = "dashed", color = "white", size = 1/2, alpha = 1/2) +
	geom_hline(aes(yintercept = age.prop * Ne), linetype = "dashed", color = "white", size = 1/2, alpha = 1/2) +
	geom_hline(aes(yintercept = age.auto * Ne), linetype = "dashed", color = "white", size = 1/2, alpha = 1/2) +
	#geom_hline(yintercept = exc$q25 * Ne, linetype = "dotted", color = "red") +
	#geom_hline(yintercept = exc$q75 * Ne, linetype = "dotted", color = "blue") +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ','), minor_breaks = c(seq(1,10,by=1), seq(10,100,by=10), seq(100,1000,by=100), seq(1000,10000,by=1000), seq(10000,100000,by=10000))) +
	scale_colour_manual(values = c("Concordant" = "forestgreen", "Discordant" = "sienna")) +
	coord_cartesian(ylim = c(0.9, 1.1e5), xlim = c(-0.01, 2.01), expand = F) +
	theme_few() +
	theme(panel.background = element_rect(fill = "grey90"),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				axis.title.x = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				legend.title = element_blank(),
				legend.justification=c(1,0), legend.position=c(0.99,1-0.99),
				legend.margin = margin(0,2,1.5,1, "mm"),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2),
				panel.grid.major.y = element_line(colour = "grey97", size = 2/3),
				panel.grid.minor.y = element_line(colour = "grey95", size = 1/2)) +
	ylab("Time (generations)")


ggsave(filename = sprintf("_plot_ccf.dev_%s_H_%s_100_%s.pdf", data, expos, near), width = 9, height = 8)




# x = rbindlist(lapply(split(x, sprintf("%d %d %d %d", x$SampleID0, x$Chr0, x$SampleID1, x$Chr1)), function(x, rec) {
# 	x$tmrca = get.tmrca(x$position[1], x$h0[1], x$h1[1], rec)
# 	x
# }, rec))

# x$tmrca = NA
# for (i in 1:nrow(x)) {
# 	if (i %% 10 == 0) cat(i, " of ", nrow(x), "\n")
# 	x$tmrca[i] = get.tmrca(x$position[i], x$h0[i], x$h1[i], rec)
# }

# x = x[order(abs(x$Shared-1), x$q50),]
# n = nrow(x)
# s = which(x$Shared == 1)
# o = which(x$Shared == 0)
# ns = length(s)
# no = length(o)
# 
# x$i = 1:n

#exc = exclude(x)







library(ggplot2)
library(ggthemes)

Ne = 2 * 10000

times = fread("~/Research/africa/dev/OutOfAfricaHapMap20.times.txt", header = T)


data = "err"
near = "NN"


marker = fread(sprintf("./%sH.marker.txt", data), header = T)


raw = fread(sprintf("new_%s_H_100_%s_raw.age.sites.txt", data, near), header = T)
adj = fread(sprintf("new_%s_H_100_%s_adj.age.sites.txt", data, near), header = T)

pair = fread(sprintf("new_%s_H_100_%s_raw.age.pairs.txt", data, near), header = T)

pair = cbind(pair, times[pair$MarkerID + 1, ])

pair = split(pair, list(pair$MarkerID))


for (tag in names(pair)) {
	x = pair[[ tag ]]


#x = pair[[ as.character(x$MarkerID[1]) ]]
#x = pair[[ sample(1:length(pair), 1) ]]


m = marker[which(marker$MarkerID == x$MarkerID[1]),]
title = sprintf("Simulated data:     Genetic map chr. = %d     Position = %s     Allele count = %d", m$Chromosome, m$Position, m$AlleleCount1)


y = rbind(cbind(raw[ which(raw$MarkerID == x$MarkerID[1]), ], fltr = "Raw"), 
					cbind(adj[ which(adj$MarkerID == x$MarkerID[1]), ], fltr = "Adjusted"))

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

x$concord = NA
x$concord[ which(x$Shared == 1) ] = "Concordant pairs"
x$concord[ which(x$Shared == 0) ] = "Discordant pairs"

x$Clock[ which(x$Clock == "M") ] = "(a) Mutation clock"
x$Clock[ which(x$Clock == "R") ] = "(b) Recombination clock"
x$Clock[ which(x$Clock == "C") ] = "(c) Combined clock"

y$Clock[ which(y$Clock == "M") ] = "(a) Mutation clock"
y$Clock[ which(y$Clock == "R") ] = "(b) Recombination clock"
y$Clock[ which(y$Clock == "C") ] = "(c) Combined clock"

ggplot(x) +
	facet_grid(.~Clock) +
	geom_rect(data = data.table(xmin = 0.95, xmax = 1.05, ymin = 1e-8, ymax = 1e8), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", alpha = 1/2) +
	geom_rect(data = data.table(xmin = -10, xmax = 10, ymin = x$node.time[1], ymax = x$parent.time[1]), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey50", alpha = 1/5) +
	geom_hline(yintercept = c(x$node.time[1], x$parent.time[1]), alpha = 1/5, size = 1/3, color = "grey40") +
	geom_linerange(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, color = concord), alpha = 1/2, size = 1/3, show.legend = F) +
	geom_point(aes(x=i, y=q50 * Ne, color = concord), shape = 16, alpha = 1, size = 1/2, show.legend = F) +
	geom_rect(data = y, aes(xmin = 0.955, xmax = 1.045, ymin=PostCI025 * Ne, ymax=PostCI975 * Ne, fill = fltr), alpha = 1/2) +
	geom_segment(data = y, aes(x = 0.956, xend = 1.044, y = PostMode * Ne, yend = PostMode * Ne, linetype = fltr), alpha = 1/2, size = 1/2) +
	scale_x_continuous(breaks = c(0.45, 1, 1.55), labels = c("Concordant\npairs", "Age\nestimate", "Discordant\npairs")) +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ','), minor_breaks = c(seq(1,10,by=1), seq(10,100,by=10), seq(100,1000,by=100), seq(1000,10000,by=1000), seq(10000,100000,by=10000))) +
	scale_colour_manual(values = c("Concordant pairs" = "forestgreen", "Discordant pairs" = "chocolate")) +
	scale_fill_manual(values = c("Raw" = "goldenrod", "Adjusted" = "dodgerblue")) +
	scale_linetype_manual(values = c("Raw" = "11", "Adjusted" = "solid")) +
	coord_cartesian(ylim = c(0.8, 1.2e5), xlim = c(-0.075, 2.075), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				panel.background = element_rect(fill = "grey92"),
				panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 12),
				axis.title.x = element_blank(),
				axis.ticks.x = element_blank(),
				axis.text.x = element_text(size = 8),
				axis.title.y = element_text(margin = margin(1,-1,1,1, "mm")),
				plot.title = element_text(size = 14),
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


ggsave(filename = sprintf("_plot_ccf.%d.new_%s_H_100_%s.pdf", x$MarkerID[1], data, near), width = 297, height = 210, units = "mm")

}






















files = dir(pattern = "^expos_.+\\.txt$")

for (file in files) {
	ex = readLines(file)
	ex = sample(ex, min(c(5000, length(ex))))
	ex = ex[order(as.numeric(ex))]
	cat(ex, sep = "\n", file = sprintf("n5000_%s", file))
}



