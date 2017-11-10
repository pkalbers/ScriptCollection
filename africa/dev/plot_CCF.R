


library(data.table)
library(ggplot2)

Ne = 2 * 7300 # 10000

marker = fread("./truH.marker.txt", header = T)
times = fread("~/Research/africa/dev/OutOfAfricaHapMap20.times.txt", header = T)
mut = as.data.table(read.table("~/Research/DetectIBD/OutOfAfricaHapMap20.mutations.txt", header = T))
rec = as.data.table(read.table("~/Research/DetectIBD/OutOfAfricaHapMap20.records.txt", header = T))

rec$c0 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\1", rec$children))
rec$c1 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\2", rec$children))

rec$children = NULL
rec$population = NULL


age = fread("test_err_fixed.age.sites.HMM.MutRec.HardBreaks.txt", header = T)

d = fread("test_err_fixed.age.pairs.HMM.MutRec.HardBreaks.txt", header = T)
d = cbind(d, times[d$MarkerID + 1, ])

d$h0 = d$SampleID0 * 2
d$h1 = d$SampleID1 * 2
x = which(d$Chr0 == "P")
d$h0[x] = d$h0[x] + 1
x = which(d$Chr1 == "P")
d$h1[x] = d$h1[x] + 1



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


d$tmrca = NA


d = split(d, d$MarkerID)


x=d$`667757`

x = d[[ sample(1:length(d), 1) ]]



for (i in 1:nrow(x)) {
	if (i %% 10 == 0) cat(i, "of", nrow(x), "\n")
	x$tmrca[i] = get.tmrca(x$position[i], x$h0[i], x$h1[i], rec)
}


x.age = age$PostMode[ which(age$MarkerID == x$MarkerID[1]) ]

x$mean = (x$Shape) / (x$Rate)
x$mode = (x$Shape-1) / x$Rate

#x$tmrca = sapply(split(x, 1:nrow(x)), function(x, r) { get.tmrca(x$position[1], x$h0[1]+1, x$h1[1]+1, r) }, rec )


x = x[order(abs(x$Shared-1), x$q50),]
n = nrow(x)
s = which(x$Shared == 1)
o = which(x$Shared == 0)
ns = length(s)
no = length(o)

x$i = 1:n

exc = exclude(x)

ggplot(x) +
	geom_pointrange(aes(x=i, y=q50 * Ne, ymin=q25 * Ne, ymax=q75 * Ne, color = factor(Shared)), alpha = 0.75) +
	geom_point(aes(x=i, y = tmrca)) +
	geom_point(aes(x=i, y = mean * Ne), colour = "black", alpha = 1/2, shape = 21) +
	geom_hline(yintercept = c(x$node.time[1], x$parent.time[1]), alpha = 0.5) +
	geom_hline(yintercept = exc$q50 * Ne, linetype = "dashed") +
	geom_hline(yintercept = x.age * Ne, linetype = "dotted") +
	#geom_hline(yintercept = exc$q25 * Ne, linetype = "dotted", color = "red") +
	#geom_hline(yintercept = exc$q75 * Ne, linetype = "dotted", color = "blue") +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ',')) +
	coord_cartesian(ylim = c(0.9e-4, 4.1) * Ne)




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


