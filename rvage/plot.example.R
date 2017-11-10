


library(ggplot2)
library(ggthemes)
library(data.table)


M = read.table("OutOfAfricaHapMap20.mutations.txt", header=T)
R = read.table("OutOfAfricaHapMap20.records.txt", header=T)

read = function(file, M, R) {
	
	age = fread(file, header = T, stringsAsFactors = F)
	
	age = as.data.frame(age)[, 1:12]
	
	age$time = NA
	
	for (i in 1:nrow(age)) {
		node = M$node[ age$MarkerID[i] + 1 ]
		
		x = which(R$node == node)
		if (length(x) == 0) next
		
		x = unique(R$time[ x ])
		
		age$time[i] = x[1]
	}
	
	if (any(is.na(age$time))) return(age[which(!is.na(age$time)), ])
	
	return(age)
}


files = "plot.age.sites.FGT.MutRec.txt"

age = NULL

for (file in files) {
	cat(file, "\n")
	tmp = read(file, M, R)
	
	tmp$Method = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.txt", "\\1", file)
	
	if (!grepl("HardBreaks", file)) {
		tmp$Type = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.txt", "\\2", file)
	} else {
		tmp$Method = sprintf("%s (HardBrks)", tmp$Method)
		tmp$Type = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.(.+)\\.txt", "\\2", file)
	}
	
	age = rbind(age, as.data.table(tmp))
}



d = age

del = which(is.na(d$PostMode))
if (length(del) > 0)
	d = d[-del, ]

del = which(is.na(d$time))
if (length(del) > 0)
	d = d[-del, ]

d$est = d$PostMode * 0.5 * 10000
d$frq = (d$Fk / 5000) * 100


sites = d
sdist = read.table("plot.age.sites.FGT.MutRec.distr.txt", header = T, stringsAsFactors = F)
pairs = read.table("plot.age.pairs.FGT.MutRec.txt", header = T, stringsAsFactors = F)
pdist = read.table("plot.age.pairs.FGT.MutRec.distr.txt", header = T, stringsAsFactors = F)

times = as.numeric(sub("^t(.+)$", "\\1", names(sdist)[-1]))
times = times * 0.5 * 10000
names(times) = names(sdist)[-1]


i = 2

mid = sites$MarkerID[i]

S = sites[which(sites$MarkerID == mid), ]
SD = sdist[which(sdist$MarkerID == mid), ]
P = pairs[which(pairs$MarkerID == mid), ]
PD = pdist[which(pdist$MarkerID == mid), ]


p = PD[, -1]

p$pair = sprintf("%d %d", P$SampleID0, P$SampleID1)
p$share = P$Shared

p = p[rev(1:nrow(p)), ]

p = melt(p, id.vars = c("pair", "share"))
names(p) = c("pair", "Shared", "Time", "CCF")
p$Shared = factor(p$Shared)
p$Time = times[p$Time]


s = SD[, -1]

s = melt(s)
names(s) = c("Time", "CCF")
s$Time = times[s$Time]



gg = ggplot(p) +
	geom_line(aes(x = Time, y = CCF, colour = Shared, group = pair), size = 0.5, alpha = 0.5) +
	geom_line(data = p[which(p$Shared == "1"), ], aes(x = Time, y = CCF, colour = Shared, group = pair), size = 0.5, alpha = 0.5) +
	geom_line(data = s, aes(x = Time, y = CCF), size = 1) +
	geom_vline(xintercept = d$time[i], size = 1.2, alpha = 1/3) +
	geom_vline(xintercept = d$time[i], size = 0.6, linetype = "dashed") +
	scale_color_manual(values = c("0"="dodgerblue", "1"="orangered")) +
	scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = (c("1", "10", "100", "1000", "10000", "100000"))) +
	coord_cartesian(xlim = c(1, 25000)) +
	theme_few() +
	theme(aspect.ratio = 1/4,
				legend.position="top")

ggsave(gg, filename = "_plot.example.pdf", width = 12, height = 6)
ggsave(gg, filename = "_plot.example.png", width = 12, height = 6)


#############











ggplot(d) +
	facet_grid(Method~Type) +
	stat_summary(aes(x=frq, y=time), fun.data=mean_se, alpha = 0.5, colour = "blue") +
	stat_summary(aes(x=frq, y=est), fun.data=mean_se, alpha = 0.75) +
	#geom_smooth(aes(x=frq, y=time), method='glm', linetype="dashed") +
	#geom_smooth(aes(x=frq, y=est), method='glm') +
	#geom_smooth(aes(x=frq, y=time), linetype="dashed") +
	#geom_smooth(aes(x=frq, y=est)) +
	coord_cartesian(ylim = c(1, 50000)) +
	theme(aspect.ratio = 1) +
	scale_y_log10() +
	xlab("Allele frequency (%)") + ylab("Mean age (generations)")



ggplot(d) +
	facet_grid(.~Type) +
	stat_summary(aes(x=frq, y=time), fun.data=mean_se, alpha = 0.25) +
	stat_summary(aes(x=frq, y=est, colour=Method), fun.data=mean_se, alpha = 0.75) +
	coord_cartesian(ylim = c(1, 50000)) +
	theme(aspect.ratio = 1) +
	scale_y_log10() +
	xlab("Allele frequency (%)") + ylab("Mean age (generations)")



ggplot(d) +
	facet_grid(Method~Type) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, est, colour = Fk)) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.5, 20000), ylim = c(0.5, 20000)) +
	theme(aspect.ratio = 1)  +
	xlab("True age") + ylab("Posterior mean age")



s = by(d, list(d$Method, d$Type, d$Fk), function(x) { if (length(x) == 0) return; cor(x$time, x$est, method = "s") } )
s =array(s, dim(s), dimnames(s))
s = melt(s)

ggplot(s) +
	facet_grid(Var1~Var2) +
	geom_point(aes(factor(Var3), value)) +
	coord_cartesian(ylim=c(0, 1)) +
	theme(aspect.ratio = 0.5) +
	xlab("fk") + ylab("Spearman rank correlation")


x = by(d, list(d$Method, d$Type), function(x) cor(x$time, x$est, method = "s") )
array(x, dim(x), dimnames(x))

x = by(d, list(d$Method, d$Type), function(x) median(abs(x$time - x$est)) )
array(x, dim(x), dimnames(x))






M = read.table("vanilla.mutations.txt", header=T)
R = read.table("vanilla.records.txt", header=T)


R$a = as.numeric(sub("^(.+),(.+)$", "\\1", R$children))
R$b = as.numeric(sub("^(.+),(.+)$", "\\2", R$children))

markers = read.table("vanilla.marker.txt", header = T, stringsAsFactors = F)

pos = markers$Position

file = "noprob.age.pairs.FGT.Rec.txt"

seg = read.table(file, header = T, stringsAsFactors = F)
seg = seg[, 1:12]
seg$trueLHS = NA
seg$trueRHS = NA

for (i in 1:nrow(seg)) {
	node = M$node[ seg$MarkerID[i] + 1 ]
	chld = R[ which(R$node == node), ]
	
	a = min(c(seg$SampleID0[i], seg$SampleID1[i]))
	b = max(c(seg$SampleID0[i], seg$SampleID1[i]))
	
	k = which(chld$a == a & chld$b == b)
	
	if (length(k) != 1) stop("???")
	
	l = unique(R$left[ which(R$node == node) ])
	r = unique(R$right[ which(R$node == node) ])
	seg$trueLHS[i] = l[1]
	seg$trueRHS[i] = r[1]
}

seg$detLHS = pos[ seg$SegmentLHS + 1 ]
seg$detRHS = pos[ seg$SegmentRHS + 1 ]



read = function(file, M, R, Ne) {
	
	age = fread(file, header = T, stringsAsFactors = F)
	
	age = as.data.frame(age)[, 1:12]
	
	age$time = NA
	
	for (i in 1:nrow(age)) {
		node = M$node[ age$MarkerID[i] + 1 ]
		
		x = which(R$node == node)
		if (length(x) == 0) next
		
		x = unique(R$time[ x ])
		#l = unique(R$left[ which(R$node == node) ])
		#r = unique(R$right[ which(R$node == node) ])
		
		age$time[i] = x[1] # / (Ne * 4)
		#age$LHS[i] = l[1]
		#age$RHS[i] = r[1]
	}
	
	if (any(is.na(age$time))) return(age[which(!is.na(age$time)), ])
	
	age
}



x=read("t1000eP.age.sites.DHG.Rec.txt", M, R, 10000)

x$PostMean = x$PostMean * 1 * 7300

ggplot(x) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, PostMean, colour = Fk)) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.5, 20000), ylim = c(0.5, 20000)) +
	theme(aspect.ratio = 1) 

cor(x$time, x$PostMean, method = "s")
cor(x$time, x$PostMode, method = "s")
median(abs(x$time - x$PostMean))

plot(sapply(split(x, x$Fk), function(x) cor(x$time, x$PostMean, method = "s")))
plot(sapply(split(x, x$Fk), function(x) median(abs(x$time - x$PostMean))))
plot(sapply(split(x, x$Fk), function(x) median(x$PostMean)))
points(sapply(split(x, x$Fk), function(x) median(x$time)), col="red")



a=read("x.age.sites.FGT.Mut.txt", M, R, 10000)
b=read("x.age.sites.FGT.Rec.txt", M, R, 10000)
c=read("x.age.sites.FGT.MutRec.txt", M, R, 10000)
z = rbind(cbind(a, Clock = "Mut"),
					cbind(b, Clock = "Rec"),
					cbind(c, Clock = "MutRec"))

ggplot(z) +
	facet_grid(.~Clock) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, PostMean, colour = Fk)) +
	scale_x_log10() +
	scale_y_log10() +
	theme(aspect.ratio = 1)


cor(a$time, a$PostMean, method = "s")
cor(b$time, b$PostMean, method = "s")
cor(c$time, c$PostMean, method = "s")



files = dir(pattern = "^error\\.H\\.age\\.sites\\..+")
files = dir(pattern = "^truth\\.P\\.age\\.sites\\..+")

files = dir(pattern = "^WE.+\\.age\\.sites\\..+")

files = dir(pattern = "^t1000eP\\.age\\.sites\\.HMM.+")
#files = files[(grep("HardBreaks", files))]

age = NULL

for (file in files) {
	cat(file, "\n")
	tmp = read(file, M, R, 10000)
	
	tmp$Method = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.txt", "\\1", file)
	
	if (!grepl("HardBreaks", file)) {
		tmp$Type = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.txt", "\\2", file)
	} else {
		tmp$Method = sprintf("%s (HardBrks)", tmp$Method)
		tmp$Type = sub(".+\\.sites\\.([A-Z]+)\\.(.+)\\.(.+)\\.txt", "\\2", file)
	}
	
	age = rbind(age, as.data.table(tmp))
}


d = age

del = which(is.na(d$PostMean))

if (length(del) > 0)
	d = d[-del, ]

d$PostMean = d$PostMean * 2 * 7300
d$frq = (d$Fk / 10000) * 100

ggplot(d) +
	facet_grid(Method~Type) +
	stat_smooth(aes(x=frq, y=time), linetype="dashed") +
	stat_smooth(aes(x=frq, y=PostMean)) +
	coord_cartesian(ylim = c(1, 10000)) +
	theme(aspect.ratio = 1) +
	scale_y_log10() +
	xlab("Allele frequency (%)") + ylab("Mean age (generations)")


ggplot(d) +
	facet_grid(Method~Type) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, PostMean, colour = Fk)) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.5, 20000), ylim = c(0.5, 20000)) +
	theme(aspect.ratio = 1)  +
	xlab("True age") + ylab("Posterior mean age")


m = by(d, list(d$Method, d$Type), function(x) { if (length(x) == 0) return; median(abs(x$time - x$PostMean)) } )
m =array(m, dim(m), dimnames(m))
m = melt(m)
names(m) = c("Method", "Type", "value")

ggplot(d) +
	facet_grid(Method~Type) +
	geom_boxplot(aes(factor(Fk), abs(time - PostMean)), outlier.colour="white", outlier.size=0) +
	#geom_hline(data = m, aes(yintercept = value), colour = "blue", alpha = 0.5) +
	coord_cartesian(ylim = c(0, 2000)) +
	theme(aspect.ratio = 1) +
	xlab("fk") + ylab("Abs. difference")



s = by(d, list(d$Method, d$Type, d$Fk), function(x) { if (length(x) == 0) return; cor(x$time, x$PostMean, method = "s") } )
s =array(s, dim(s), dimnames(s))
s = melt(s)

S = by(d, list(d$Method, d$Type), function(x) { if (length(x) == 0) return; cor(x$time, x$PostMean, method = "s") } )
S =array(S, dim(S), dimnames(S))
S = melt(S)

ggplot(s) +
	facet_grid(Var1~Var2) +
	geom_point(aes(factor(Var3), value)) +
	geom_hline(data = S, aes(yintercept = value), colour = "blue", alpha = 0.5) +
	coord_cartesian(ylim=c(0, 1)) +
	theme(aspect.ratio = 1) +
	xlab("fk") + ylab("Spearman rank correlation")


x = by(d, list(d$Method, d$Type), function(x) cor(x$time, x$PostMean, method = "s") )
array(x, dim(x), dimnames(x))

x = by(d, list(d$Method, d$Type), function(x) median(abs(x$time - x$PostMean)) )
array(x, dim(x), dimnames(x))





x = by(d, list(d$Method, d$Type), function(x) cor(x$time, x$PostMean, method = "p")^2 )
array(x, dim(x), dimnames(x))
x = by(d, list(d$Method, d$Type), function(x) cor(x$time, x$PostMean, method = "s") )
array(x, dim(x), dimnames(x))

x = by(d, list(d$Method, d$Type), function(x) (sqrt(mean(((log(x$time) - log(x$PostMean))^2)))) )
array(x, dim(x), dimnames(x))

x = by(d, list(d$Method, d$Type), function(x) median(abs(x$time - x$PostMean)) )
array(x, dim(x), dimnames(x))


d$delta = abs(d$time - d$PostMean)

ggplot(d) +
	#facet_grid(Method~Type) +
	facet_grid(.~Data) +
	geom_boxplot(aes(Method, delta)) +
	coord_cartesian(ylim = c(0, 5000))





mr = "./truth.H.age.sites.FGT.MutRec.txt"
we = "./WE.truth.H.age.sites.FGT.MutRec.txt"

mr = read(mr, M, R, 7300)
we = read(we, M, R, 7300)

mr$Theta = "4 Ne µ"
we$Theta = "Watterson"

d = rbind(mr, we)


x = by(d, d$Theta, function(x) cor(x$time, x$PostMean, method = "p")^2 )
array(x, dim(x), dimnames(x))

x = by(d, d$Theta, function(x) sqrt(mean(((x$time - x$PostMean)^2))) )
array(x, dim(x), dimnames(x))

x = by(d, d$Theta, function(x) median(abs(x$time - x$PostMean)) )
array(x, dim(x), dimnames(x))


ggplot(d) +
	facet_grid(.~Theta) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, PostMean, colour = Fk)) +
	scale_x_log10() +
	scale_y_log10()




function (hmm, observation) 
{
	hmm$transProbs[is.na(hmm$transProbs)] = 0
	hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
	
	nObservations = length(observation)
	nStates = length(hmm$States)
	
	v = array(NA, c(nStates, nObservations))
	dimnames(v) = list(states = hmm$States, index = 1:nObservations)
	
	for (state in hmm$States) {
		v[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, observation[1]])
	}
	for (k in 2:nObservations) {
		for (state in hmm$States) {
			maxi = NULL
			for (previousState in hmm$States) {
				temp = v[previousState, k - 1] + log(hmm$transProbs[previousState, state])
				maxi = max(maxi, temp)
			}
			v[state, k] = log(hmm$emissionProbs[state, observation[k]]) + maxi
		}
	}
	viterbiPath = rep(NA, nObservations)
	for (state in hmm$States) {
		if (max(v[, nObservations]) == v[state, nObservations]) {
			viterbiPath[nObservations] = state
			break
		}
	}
	for (k in (nObservations - 1):1) {
		for (state in hmm$States) {
			if (max(v[, k] + log(hmm$transProbs[, viterbiPath[k + 1]])) == v[state, k] + log(hmm$transProbs[state, viterbiPath[k + 1]])) {
				viterbiPath[k] = state
				break
			}
		}
	}
	return(viterbiPath)
}


