


library(data.table)
library(ggplot2)
library(ggthemes)



M = read.table("OutOfAfricaHapMap20.mutations.txt", header=T)
R = read.table("OutOfAfricaHapMap20.records.txt", header=T)

read.age = function(file, M, R) {
	
	age = fread(file, header = T, stringsAsFactors = F)
	
	age$node = M$node[ age$MarkerID + 1 ]
	
	age$time = match(age$node, R$node)
	
	del = which(is.na(age$time))
	if (length(del) > 0)
		age = age[-del, ]
	
	age$time = R$time[ age$time ]
	
	# for (i in 1:nrow(age)) {
	# 	x = which(R$node == age$node[i])
	# 	if (length(x) == 0) next
	# 	x = unique(R$time[ x ])
	# 	age$time[i] = x[1]
	# }
	# if (any(is.na(age$time))) return(age[which(!is.na(age$time)), ])
	
	return(age)
}




markers = fread("truH.marker.txt", header = T, stringsAsFactors = F)

files = dir(pattern = "^x100\\.tru.+\\.age\\.sites\\..+")


d = NULL

for (file in files) {
	
	method = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\1", file)
	clock  = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\2", file)
	
	phased = grepl(".+P\\.age\\..+", file)
	
	if (phased)
		method = sprintf("%s phased", method)
	
	#if (clock == "MutRec") clock = " MutRec"
	
	cat(".")
	
	age = read.age(file, M, R)
	
	age$method = method
	age$clock = clock
	
	d = rbind(d, age)
}


del = which(d$N_Others < d$N_Shared)
if (length(del) > 0)
	d = d[-del, ]


m = match(d$MarkerID, markers$MarkerID)
d$a = markers$AlleleCount1[m]
d$g = markers$GenotypeCount1[m]
del = which(d$a != d$g)
if (length(del) > 0)
	d = d[-del, ]


#d = d[sample(1:nrow(d), 50000), ]


d$frq = (d$Fk / 5000) * 100

d$est = NA

x = (d$clock == "MutRec")
d$est[x] = d$PostMean[x] * 7300 * (2/100) * 2
x = (d$clock == "Mut")
d$est[x] = d$PostMean[x] * 7300 * (2/100)
x = (d$clock == "Rec")
d$est[x] = d$PostMean[x] * 7300 * 2

x = (d$method == "FGT" | d$method == "FGT phased")
d$est[x] = d$PostMean[x] * 7300 * (1/4)


x = by(d, list(d$clock, d$method), function(x) median(abs(x$time - x$est)) )
array(x, dim(x), dimnames(x))

x = by(d, list(d$clock, d$method), function(x) cor(x$time, x$est, method = "s") )
array(x, dim(x), dimnames(x))



gg = ggplot(d) +
	facet_grid(clock~method) +
	geom_point(aes(x=frq, y=est), colour = "grey50", alpha = 0.05, shape=15, size=1) +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=2.5, colour = "white") +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=1, colour = "grey10") +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=1.5, colour = "blue") +
	coord_cartesian(ylim = c(1, 2000)) +
	scale_y_log10(breaks=c(0, 1, 10, 100, 1000, 10000, 100000)) +
	scale_x_continuous(breaks = (0:5)/10) +
	theme_few() +
	theme(aspect.ratio=2.1/3) +
	xlab("Allele frequency (%)") + ylab("Inferred age (generations)")

ggsave(gg, filename = "_plot.ooa.tru.png", width = 12, height = 10)
ggsave(gg, filename = "_plot.ooa.err.png", width = 12, height = 10)



ggplot(d) +
	facet_grid(clock~method) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, est, colour = Fk), size = 0.5) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.05, 50000), ylim = c(0.05, 50000)) +
	theme(aspect.ratio = 1)  +
	xlab("True age") + ylab("Posterior mode age")





files = dir(pattern = "^t1000\\.truH\\.age\\.pairs\\..+")

p = lapply(files, function(file) {
	x = fread(file, header = T, stringsAsFactors = F)
	x$method = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\1", file)
	x$clock  = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\2", file)
	x
})
p = rbindlist(p)

p = p[which(p$Shared == 1),]
p = p[-(which(p$SegmentLHS == 0)),]
p = p[-(which(p$SegmentRHS == max(p$SegmentRHS))),]

pos = markers$Position
names(pos) = as.character(markers$MarkerID)

p$pos_SegmentLHS = pos[ as.character(p$SegmentLHS) ]
p$pos_SegmentRHS = pos[ as.character(p$SegmentRHS) ]

p$len = (p$pos_SegmentRHS - p$pos_SegmentLHS) + 1

p$str = sprintf("%d %d %d", p$MarkerID + 1, p$SampleID0 + 1, p$SampleID1 + 1)



load("~/Research/DetectIBD/result.truth.local.RData")

truth$str = sprintf("%d %d %d", truth$index, truth$g0, truth$g1)


p = split(p, list(p$method, p$clock))

for (tag in names(p)) {
	sub = p[[tag]]
	
	i = intersect(sub$str, truth$str)
	
	a = match(i, sub$str)
	b = match(i, truth$str)
	
	sub = sub[a, ]
	t = truth[b, ]
	
	sub$len_true = (t$rhs.position - t$lhs.position) + 1
	
	p[[tag]] = sub
}

p = rbindlist(p)



ggplot(p[sample(1:nrow(p), 100000),]) +
	facet_grid(method~clock) +
	geom_point(aes(len_true, len)) +
	scale_x_log10() +
	scale_y_log10()
	stat_summary(aes(x=Fk, y=time), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	geom_boxplot(aes(factor(Fk), len))
	






