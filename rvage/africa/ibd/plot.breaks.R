#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)
library(colorspace)


load("~/Research/DetectIBD/GenErr_1000G/result.match.OutOfAfricaHapMap20.GenErr_1000G.RData")
TR = as.data.table(match)


del = which(TR$true.lhs.idx >= TR$true.rhs.idx); if (length(del) > 0) TR = TR[-del,]

del = which(TR$true.lhs.idx >= TR$index); if (length(del) > 0) TR = TR[-del,]
del = which(TR$true.rhs.idx <= TR$index); if (length(del) > 0) TR = TR[-del,]


M = fread("../truH.marker.txt")
POS = M$Position



load("~/Research/DetectIBD/results.hmm_tru.RData")
dn = data; rm(data)
dt = fread("all.truG.ibd.HMM.txt", header = T, stringsAsFactors = F)  ##  11315116
de = fread("all.errG.ibd.HMM.txt", header = T, stringsAsFactors = F)  ##  11035589


del = which(dn$LHS >= dn$RHS); if (length(del) > 0) dn = dn[-del,]
del = which(dt$LHS >= dt$RHS); if (length(del) > 0) dt = dt[-del,]
del = which(de$LHS >= de$RHS); if (length(del) > 0) de = de[-del,]

del = which(dn$LHS >= dn$MarkerID); if (length(del) > 0) dn = dn[-del,]
del = which(dn$RHS <= dn$MarkerID); if (length(del) > 0) dn = dn[-del,]

del = which(dt$LHS >= dt$MarkerID); if (length(del) > 0) dt = dt[-del,]
del = which(dt$RHS <= dt$MarkerID); if (length(del) > 0) dt = dt[-del,]

del = which(de$LHS >= de$MarkerID); if (length(del) > 0) de = de[-del,]
del = which(de$RHS <= de$MarkerID); if (length(del) > 0) de = de[-del,]


dn$key = sprintf("%d %d %d", dn$MarkerID+1, dn$SampleID0+1, dn$SampleID1+1)
dt$key = sprintf("%d %d %d", dt$MarkerID+1, dt$SampleID0+1, dt$SampleID1+1)
de$key = sprintf("%d %d %d", de$MarkerID+1, de$SampleID0+1, de$SampleID1+1)
TR$key = sprintf("%d %d %d", TR$index, TR$g0, TR$g1)

any(duplicated(dn$key))
any(duplicated(dt$key))
any(duplicated(de$key))
any(duplicated(TR$key))

xn = intersect(TR$key, dn$key)   ## 
xt = intersect(TR$key, dt$key)   ## 
xe = intersect(TR$key, de$key)   ## 
key = intersect(intersect(xt, xe), xn)  ## 

dn = as.data.table(as.data.frame(dn)[which(dn$key %in% key),])
dt = as.data.table(as.data.frame(dt)[which(dt$key %in% key),])
de = as.data.table(as.data.frame(de)[which(de$key %in% key),])
TR = as.data.table(as.data.frame(TR)[which(TR$key %in% key),])

dn = dt[order(dn$key), ]
dt = dt[order(dt$key), ]
de = de[order(de$key), ]
TR = TR[order(TR$key), ]

identical(dt$key, de$key)
identical(dt$key, TR$key)
identical(dn$key, TR$key)



match = data.table(index = TR$index, pos = TR$pos, fk = TR$fk, g0=TR$g0, g1=TR$g1, h0=TR$h0, h1=TR$h1, time=TR$true.time, 
									 n.lhs = dn$LHS+1,
									 n.rhs = dn$RHS+1,
									 t.lhs = dt$LHS+1,
									 t.rhs = dt$RHS+1,
									 e.lhs = de$LHS+1,
									 e.rhs = de$RHS+1,
									 true.lhs = TR$true.lhs.idx,
									 true.rhs = TR$true.rhs.idx,
									 pos.n.lhs = M$Position[dn$LHS+1],
									 pos.n.rhs = M$Position[dn$RHS+1],
									 pos.t.lhs = M$Position[dt$LHS+1],
									 pos.t.rhs = M$Position[dt$RHS+1],
									 pos.e.lhs = M$Position[de$LHS+1],
									 pos.e.rhs = M$Position[de$RHS+1],
									 pos.true.lhs = M$Position[TR$true.lhs.idx],
									 pos.true.rhs = M$Position[TR$true.rhs.idx])

save(match, file="match.neu_tru_err.RData")
load("match.neu_tru_err.RData")

save(match, file="match.true_err.RData")
load("match.true_err.RData")



### stats


# number of sites per fk
x = by(match$index, match$fk, function(x) length(unique(x)))
x = array(x, dim(x), dimnames(x))


# average per pair
k = sprintf("%d %d", match$g0, match$g1)
k = table(k)
mean(k); se(k)




wl = POS[2]
wr = POS[length(POS)-1]
nr = nrow(match)


match = match[order(match$fk),]

# unique
key0 = sprintf("%d %d %d %d", match$n.lhs, match$n.rhs, match$g0, match$g1);  length(unique(key0)) / nr * 100
key1 = sprintf("%d %d %d %d", match$t.lhs, match$t.rhs, match$g0, match$g1);  length(unique(key1)) / nr * 100
key2 = sprintf("%d %d %d %d", match$e.lhs, match$e.rhs, match$g0, match$g1);  length(unique(key2)) / nr * 100
key4 = sprintf("%d %d %d %d", match$true.lhs, match$true.rhs, match$g0, match$g1);  length(unique(key4)) / nr * 100

del1 = duplicated(key1)
del2 = duplicated(key2)
del4 = duplicated(key4)

del = unique(c(which(del1), which(del2), which(del4)))
((length(del) / nr)) * 100
match = match[-del, ]
nrow(match) / nr * 100
nr = nrow(match)


# tmp = split(match, match$fk)
# tmp = lapply(tmp, function(match) {
# 	key1 = sprintf("%d %d %d %d", match$t.lhs, match$t.rhs, match$g0, match$g1);
# 	key2 = sprintf("%d %d %d %d", match$e.lhs, match$e.rhs, match$g0, match$g1);
# 	key4 = sprintf("%d %d %d %d", match$true.lhs, match$true.rhs, match$g0, match$g1);
# 	
# 	del1 = duplicated(key1)
# 	del2 = duplicated(key2)
# 	del4 = duplicated(key4)
# 	
# 	del = unique(c(which(del1), which(del2), which(del4)))
# 	match[-del, ]
# })
# tmp = rbindlist(tmp)
# match=tmp


# average per pair
k = sprintf("%d %d", match$g0, match$g1)
k = table(k)
mean(k); se(k)


# bounds

length(which((match$h.lhs <= wl & match$h.rhs < wr) | (match$h.lhs > wl & match$h.rhs >= wr))) / nr * 100
length(which(match$h.lhs <= wl & match$h.rhs >= wr)) / nr * 100

length(which((match$p.lhs <= wl & match$p.rhs < wr) | (match$p.lhs > wl & match$p.rhs >= wr))) / nr * 100
length(which(match$p.lhs <= wl & match$p.rhs >= wr)) / nr * 100

length(which((match$g.lhs <= wl & match$g.rhs < wr) | (match$g.lhs > wl & match$g.rhs >= wr))) / nr * 100
length(which(match$g.lhs <= wl & match$g.rhs >= wr)) / nr * 100

length(which((match$true.lhs.idx == 1 & match$true.rhs.idx != length(POS)) | (match$true.lhs.idx != 1 & match$true.rhs.idx == length(POS)))) / nr * 100
length(which(match$true.lhs.idx == 1 & match$true.rhs.idx == length(POS))) / nr * 100


match$t.wall = (match$t.lhs <= 2 | match$t.rhs >= nrow(M)-2)  ;  t=table(match$t.wall); t/sum(t)*100
match$e.wall = (match$e.lhs <= 2 | match$e.rhs >= nrow(M)-2)  ;  t=table(match$e.wall); t/sum(t)*100
match$true.wall = (match$true.lhs <= 2 | match$true.rhs >= nrow(M)-2)  ;  t=table(match$true.wall); t/sum(t)*100



# del = which(match$g.wall | match$h.wall | match$p.wall | match$true.wall) # | match$hmm.wall)
# length(del) / nr * 100
# if (length(del) > 0) {
# 	match = match[-del, ]
# 	nr = nrow(match)
# }



d = rbind(data.table(side="LHS", fk = match$fk, type = "(i) HMM, before error, neutral model",  det = match$pos - match$pos.n.lhs + 1,   tru = match$pos - match$pos.true.lhs + 1),
					data.table(side="LHS", fk = match$fk, type = "(ii) HMM, before error, empirical model",  det = match$pos - match$pos.t.lhs + 1,   tru = match$pos - match$pos.true.lhs + 1),
					data.table(side="LHS", fk = match$fk, type = "(iii) HMM, after error, empirical model",   det = match$pos - match$pos.e.lhs + 1,   tru = match$pos - match$pos.true.lhs + 1),
					data.table(side="RHS", fk = match$fk, type = "(i) HMM, before error, neutral model",  det = match$pos.n.rhs - match$pos + 1,   tru = match$pos.true.rhs - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(ii) HMM, before error, empirical model",  det = match$pos.t.rhs - match$pos + 1,   tru = match$pos.true.rhs - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "(iii) HMM, after error, empirical model",   det = match$pos.e.rhs - match$pos + 1,   tru = match$pos.true.rhs - match$pos + 1)
)

d = rbind(data.table(side="LHS", fk = match$fk, type = "HMM, before error",  det = match$pos - match$pos.t.lhs + 1,   tru = match$pos - match$pos.true.lhs + 1),
					data.table(side="LHS", fk = match$fk, type = "HMM, after error",   det = match$pos - match$pos.e.lhs + 1,   tru = match$pos - match$pos.true.lhs + 1),
					data.table(side="RHS", fk = match$fk, type = "HMM, before error",  det = match$pos.t.rhs - match$pos + 1,   tru = match$pos.true.rhs - match$pos + 1),
					data.table(side="RHS", fk = match$fk, type = "HMM, after error",   det = match$pos.e.rhs - match$pos + 1,   tru = match$pos.true.rhs - match$pos + 1)
)

table(d$type)

d$type = factor(d$type, levels = c("(i) HMM, before error, neutral model", "(ii) HMM, before error, empirical model", "(iii) HMM, after error, empirical model"), ordered = T)
d$type = factor(d$type, levels = c("HMM, before error", "HMM, after error"), ordered = T)



### STATS 
se <- function(x) sqrt(var(x)/length(x))

rmsle = function(a, b) {
	n = length(a)
	a = log10(a + 1)
	b = log10(b + 1)
	sqrt(sum((a - b)^2) / n)
}


# over/underest
ou = data.table(r = round(match$pos.e.rhs / 10), tr = round(match$pos.true.rhs / 10),
								l = round(match$pos.e.lhs / 10), tl = round(match$pos.true.lhs / 10))

r = which(ou$r > ou$tr)
l = which(ou$l < ou$tl)
(length(l)+length(r)) / (nrow(match)*2) * 100

r = which(ou$r < ou$tr)
l = which(ou$l > ou$tl)
(length(l)+length(r)) / (nrow(match)*2) * 100

r = which(ou$r == ou$tr)
l = which(ou$l == ou$tl)
(length(l)+length(r)) / (nrow(match)*2) * 100





# average per pair
k = sprintf("%d %d", match$g0, match$g1)
k = table(k)
mean(k); se(k)


x = by(abs(d$tru - d$det)/1e6, list(d$type, d$fk), mean)
array(x, dim(x), dimnames(x))

x = by(abs(d$tru - d$det)/1e6, list(d$type, d$fk), se)
array(x, dim(x), dimnames(x))

by(d, d$type, function(x) cor(x$tru, x$det, method = "p")^2)
by(d, list(d$type, d$side), function(x) cor(x$tru, x$det, method = "p")^2)
by(d, list(d$type, d$side), function(x) cor(x$tru, x$det, method = "s"))

by(d, d$type, function(x) { q = table(x$tru < x$det); (q/sum(q))*100 })
by(d, c(d$side, d$type), function(x) { q = table(x$tru < x$det); (q/sum(q))*100 })


x = by(d, list(d$type), function(x) cor(x$tru, x$det, method = "p")^2)
t(array(x, dim(x), dimnames(x)))

x = by(d, list(d$type, d$fk), function(x) { cor(x$tru, x$det, method = "p")^2 })
x = t(array(x, dim(x), dimnames(x)))

x = by(q, list(q$type, q$fk), function(x) { cor(x$tru, x$det, method = "p")^2 })
x = t(array(x, dim(x), dimnames(x)))


y = by(d, list(d$type), function(x) rmsle(x$tru, x$det))
t(array(y, dim(y), dimnames(y)))

y = by(d, list(d$type, d$fk), function(x) { rmsle(x$tru, x$det) })
y = t(array(y, dim(y), dimnames(y)))

y = by(q, list(q$type, q$fk), function(x) { rmsle(x$tru, x$det) })
y = t(array(y, dim(y), dimnames(y)))

cbind(x, y)

z = by(d, list(d$type, d$fk), function(x) nrow(x)/2)
t(array(z, dim(z), dimnames(z)))





### hist, one-sided

d$map = (d$det / d$tru)

p = d

# del = which(p$map > 10)
# if (length(del) > 0) p = p[-del,]

k = which(p$fk %in% c(2, 5, 10, 15, 20, 25))
p = p[k, ]

p = split(p, list(p$type, p$fk))
p = lapply(p, function(z) {
	x = seq(0, 10, by=0.05)
	y = ecdf(z$map)(x)
	data.table(type = z$type[1], fk = z$fk[1], x=x, y=y)
})
p = rbindlist(p)

p$fk = factor(p$fk, levels = c(2:25), ordered = T)

hist = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	geom_vline(xintercept=c(1), colour="white", size = 1.25) +
	geom_line(aes(x, y, color = fk), size = 1, alpha = 0.9) +
	geom_vline(xintercept=c(1), colour="black", size = 0.75, linetype = "22") +
	scale_y_continuous(breaks = seq(0, 1, by=0.2), minor_breaks = seq(0, 1, by=0.1)) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
	coord_cartesian(xlim = c(-0.005, 2.7505), ylim = c(-0.005, 1.005), expand = F) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey90", "purple"))(3)[2], 
																										 colorRampPalette(c("grey70", "navy"))(3)[2], 
																										 colorRampPalette(c("grey20", "red"))(3)[2], 
																										 colorRampPalette(c("grey95", "goldenrod4"))(3)[2])))(6)) +
	theme_few() +
	theme(aspect.ratio=1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
				panel.grid.minor.y = element_line(size = 0.5, colour = "grey80"),
				panel.grid.major.x = element_line(size = 0.5, colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance to true breakpoint") +
	ylab("CDF") +
	guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
hist

ggsave(filename = "___plot.break.hist.pdf", plot = hist, width = 9, height = 6.1)





### scatter


brk = seq(log10(1), log10(100e6), length.out = 151)

p = lapply(split(d, list(d$type, d$side)), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	#mat$fk  = x$fk[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)


x = (p$side == "LHS")
p$Var1[x] = p$Var1[x] * -1


tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


tm = rbind(data.table(x = c(tik*-1 -1,tik+1), y = 0, xend = c(tik*-1 -1,tik+1), yend = 5),
					 data.table(y = tik, x = 150, yend = c(tik,tik), xend = 145),
					 data.table(y = tik, x = -150, yend = c(tik,tik), xend = -145))


scatter = ggplot(data = p) +
	facet_wrap(~type, ncol = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	geom_abline(slope = c(-1, 1), intercept = c(0,0), colour="black", alpha = 0.25) +
	#scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 2.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	#scale_fill_gradient2(low = "grey", mid = "darkblue", high = "orange", na.value = "grey85", midpoint = 1.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_fill_gradientn(colours = c("white", "snow", "yellow", "gold", "orange", "orangered1", "orangered4", "firebrick4", "firebrick4"), na.value = "grey85", limits = c(0, 5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	annotate("text", label = "LHS", x = length(brk) *-0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	annotate("text", label = "RHS", x = length(brk) * 0.1, y = length(brk) * 0.95, size = 3.5, colour="grey20") +
	scale_x_continuous(expand = c(0, 0), breaks = c(rev(-1 * tik) - 1, tik + 1), labels = c(rev(tag), tag)) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 200/401, 
				panel.background = element_rect(fill = "grey60"),
				legend.background = element_rect(fill = "grey85"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				legend.title = element_blank()) +
	ylab("Detected breakpoint distance") +
	xlab("True breakpoint distance")
scatter

ggsave(filename = "___plot.break.scat.pdf", plot = scatter, width = 9, height = 6.1)





###
### length
###

del = which(match$t.wall | match$e.wall | match$true.wall)
if (length(del) > 0) match = match[-del, ]

t.phy = ((match$pos.t.rhs - match$pos.t.lhs)+1) / 1e6
e.phy = ((match$pos.e.rhs - match$pos.e.lhs)+1) / 1e6
z.phy = ((match$pos.true.rhs - match$pos.true.lhs)+1) / 1e6

t.gen = (M$GenDist[ match$t.rhs ] - M$GenDist[ match$t.lhs ]) + 1e-08
e.gen = (M$GenDist[ match$e.rhs ] - M$GenDist[ match$e.lhs ]) + 1e-08
z.gen = (M$GenDist[ match$true.rhs ] - M$GenDist[ match$true.lhs ]) + 1e-08


q = rbind(data.table(idx = match$index, fk = match$fk, type = "HMM, before error", x = t.phy, mode = "A. Physical length"),
					data.table(idx = match$index, fk = match$fk, type = "HMM, before error", x = t.gen, mode = "B. Genetic length"),
					data.table(idx = match$index, fk = match$fk, type = "HMM, after error",  x = e.phy, mode = "A. Physical length"),
					data.table(idx = match$index, fk = match$fk, type = "HMM, after error",  x = e.gen, mode = "B. Genetic length"),
					data.table(idx = match$index, fk = match$fk, type = "True IBD",          x = z.phy, mode = "A. Physical length"),
					data.table(idx = match$index, fk = match$fk, type = "True IBD",          x = z.gen, mode = "B. Genetic length"))

q$type = factor(q$type, levels = c("True IBD", "HMM, before error", "HMM, after error"), ordered = T)


# q = split(q, list(q$type, q$mode))
# cor(q$`True IBD.A. Physical length`$x, q$`(a) FGT, true haplotypes.A. Physical length`$x, method = 'p')^2
# cor(q$`True IBD.A. Physical length`$x, q$`(b) FGT, phased haplotypes.A. Physical length`$x, method = 'p')^2


x = by(q$x, list(q$type, q$mode), median)
array(x, dim(x), dimnames(x))

x = by(q$x, list(q$fk, q$type, q$mode), median)
array(x, dim(x), dimnames(x))


q = split(q, list(q$fk, q$type, q$mode))


# sub = sample(match$index, 10000)
# q = lapply(q, function(x) {
# 	k = which(x$idx %in% sub)
# 	x[k, ]
# })
# p = rbindlist(q)

q = lapply(q, function(x) {
	z = x
	w = fivenum(x$x)
	z$low = w[2]
	z$med = w[3]
	z$upp = w[4]
	z[1,]
})
q = rbindlist(q)

#p$fk = factor(p$fk, levels = (2:25), ordered = T)
#q$fk = factor(q$fk, levels = (2:25), ordered = T)

tl = c("True IBD", "HMM, before error", "HMM, after error")

q$type = factor(q$type, levels = tl, ordered = T)

ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))

len = ggplot(q) + #[sample(1:nrow(p), 10000),]) +
	facet_wrap(~mode, scales = "free", ncol = 1) +
	#geom_violin(aes(fk, x, fill = type), size = 0, alpha = 0.9) +
	#geom_pointrange(data=q, aes(x=fk, y=med, ymin=low, ymax=upp, group=type), color = "black", size = 2/3, shape='|', position=position_dodge(width = 0.9)) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=type), size = 1.5, position=position_dodge(width = 0.8)) +
	geom_point(aes(x=fk, y=med, color=type), size = 3, shape=18, position=position_dodge(width = 0.8)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), minor_breaks = mb) +
	scale_x_continuous(breaks = 2:25) +
	scale_color_manual(values = c("grey20", "dodgerblue", "orangered")) +
	#scale_fill_manual(values = c("turquoise3", "purple", "limegreen", "grey50")) +
	coord_cartesian(ylim = c(0.025, 12.5), xlim = c(1.5, 25.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "vertical",
				legend.text = element_text(size = 9),
				legend.key.height = unit(0.1, "cm"),
				legend.background = element_rect(fill = "white", colour = "grey50", size = 0.5),
				#panel.spacing.x = unit(-1, "points"), 
				strip.text = element_text(face = "bold", hjust = 0),
				#strip.text.y = element_blank(),
				#axis.ticks.y = element_blank(),
				panel.grid = element_line(colour = "grey80", size = 1/3),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") +
	ylab("Genetic length (cM)                                     Physical length (Mb)") +
	guides(color = guide_legend(override.aes = list(size = 5)))
len

ggsave(len, filename = "___plot.length.pdf", height=7, width = 9)



##   stats


x = by(q$x, list(q$type, q$mode), median)
array(x, dim(x), dimnames(x))




p = split(p, list(p$type, p$mode))

for (tag in names(p)) cat(tag, lm(x ~ fk, p[[tag]])$coefficients, "\n")



#####################



























####

d = rbind(data.table(fk = match$fk, type = "(c) DGT",                    det = (match$g.rhs - match$g.lhs)+1,     tru = (match$true.rhs.pos - match$true.lhs.pos)+1),
					data.table(fk = match$fk, type = "(a) FGT, true haplotypes",   det = (match$h.rhs - match$h.lhs)+1,     tru = (match$true.rhs.pos - match$true.lhs.pos)+1),
					data.table(fk = match$fk, type = "(b) FGT, phased haplotypes", det = (match$p.rhs - match$p.lhs)+1,     tru = (match$true.rhs.pos - match$true.lhs.pos)+1)
					#data.table(fk = match$fk, type = "Hidden Markov Model, genotype data",   det = (match$hmm.rhs - match$hmm.lhs)+1, tru = (match$true.rhs.pos - match$true.lhs.pos)+1))
)

### OR ###

d = rbind(data.table(fk = match.H$fk, type = "(a) Beagle IBD, true haplotypes",   det = match.H$beagle.rhs.position - match.H$beagle.lhs.position, tru = POS[match.H$rhs.index] - POS[match.H$lhs.index]),
					data.table(fk = match.P$fk, type = "(b) Beagle IBD, phased haplotypes", det = match.P$beagle.rhs.position - match.P$beagle.lhs.position, tru = POS[match.P$rhs.index] - POS[match.P$lhs.index]))



del = which(d$tru == 0)
if (length(del) > 0)
	d = d[-del, ]

#d$det = d$det + 1
#d$tru = d$tru + 1


### STATS 

by(d, d$type, function(x) cor(x$tru, x$det, method = "p")^2)
by(d, d$type, function(x) cor(x$tru, x$det, method = "s"))


# scatter

brk = seq(log10(1), log10(100e6), length.out = 201)

p = lapply(split(d, d$type), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)



tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


scatter = ggplot(data = p) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = 1, intercept = 0, colour="black", alpha=0.25) +
	scale_fill_gradientn(colours = c("white", "yellow", "orangered", "red"), na.value = "grey75", limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.background = element_rect(fill = "grey90"),
				axis.text = element_text(size = 8),
				legend.position = "top",
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				legend.title = element_blank()) +
	ylab("Detected segment length") +
	xlab("True segment length")

ggsave(filename = sprintf("_plot.true-inferred.length.scatter.%s.pdf", prefix), plot = scatter, width = 9, height = 9)



# boxplot

tru = lapply(split(d, list(d$type, d$fk)), function(x) {
	data.table(fk = x$fk[1], md = median(x$tru), type = x$type[1])
})
tru = rbindlist(tru)


a = sort(unique(d$fk))
b = (a / 5000) * 100
frq = b
names(frq) = a
x = ((a %% 5 == 0))
frq[ !x ] = ""


box = ggplot(d) + 
	facet_grid(.~type) + 
	geom_boxplot(aes(x = factor(fk), y = det), fill = "grey70", outlier.colour = "white", outlier.size = 0) +
	geom_point(data = tru, aes(x = factor(fk), y = md), colour = "white", size = 2.5, alpha = 0.5) +
	geom_point(data = tru, aes(x = factor(fk), y = md), colour = "blue", size = 1) +
	scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 3e7, by = 0.25e7), labels = sprintf("%.1f", seq(0, 3e7, by = 0.25e7) / 1e6)) +
	scale_x_discrete(labels = frq) +
	coord_cartesian(ylim=c(0, 1e7)) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (Mb)")







#### BEAGLE

# scatter

brk = seq(log10(1), log10(100e6), length.out = 201)

p = lapply(split(d, d$type), function(x, brk) {
	det = as.numeric(cut(log10(abs(x$det)), breaks = brk, include.lowest = T))
	tru = as.numeric(cut(log10(abs(x$tru)), breaks = brk, include.lowest = T))
	mat = matrix(0, nrow = length(brk)-1, ncol = length(brk)-1)
	for (i in 1:nrow(x)) {
		mat[ tru[i], det[i] ] = mat[ tru[i], det[i] ] + 1
	}
	mat = melt(mat, value.name = "x")
	mat$type  = x$type[1]
	mat$side  = x$side[1]
	as.data.table(mat)
}, brk)
p = rbindlist(p)



tag = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")
tik = as.numeric(cut(log10(10^(1:7)), breaks = brk, include.lowest = T))


scatter = ggplot(data = p) +
	facet_wrap(~type, nrow = 1) +
	geom_raster(aes(x = Var1, y = Var2, fill = log10(x))) +
	geom_abline(slope = 1, intercept = 0, colour="black", alpha=0.25) +
	scale_fill_gradient2(low = "darkblue", mid = "darkorange", high = "white", na.value = "grey85", midpoint = 2.5, limits = c(0, 4.5), breaks = log10(10^(0:6)), labels = sprintf("%d", 10^(0:6))) +
	scale_x_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	scale_y_continuous(expand = c(0, 0), breaks = tik, labels = tag) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold"),
				legend.position = "top",
				legend.title = element_blank()) +
	ylab("Detected segment length") +
	xlab("True segment length")

ggsave(filename = sprintf("__plot.beagle.break-true.length-scatter.%s.pdf", prefix), plot = scatter, height = 10, width = 12-4)
ggsave(filename = sprintf("__plot.beagle.break-true.length-scatter.%s.png", prefix), plot = scatter, height = 10, width = 12-4)


# boxplot

tru = lapply(split(d, list(d$type, d$fk)), function(x) {
	data.table(fk = x$fk[1], md = median(x$tru), type = x$type[1])
})
tru = rbindlist(tru)

gg = ggplot(d) + 
	facet_grid(.~type) + 
	geom_boxplot(aes(x = factor(fk), y = det), fill = "grey70", outlier.colour = "white", outlier.size = 0) +
	geom_point(data = tru, aes(x = factor(fk), y = md), colour = "white", size = 2.5, alpha = 0.5) +
	geom_point(data = tru, aes(x = factor(fk), y = md), colour = "blue", size = 1) +
	scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 3e7, by = 0.25e7), labels = sprintf("%.1f", seq(0, 3e7, by = 0.25e7) / 1e6)) +
	scale_x_discrete(labels = frq) +
	coord_cartesian(ylim=c(0, 1e7)) +
	theme_few() +
	theme(aspect.ratio = 1,
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency (%)") + ylab("Detected segment length (Mb)")

gg

ggsave(gg, filename = "__boxplot.beagle.tru.pdf", height = 10, width = 12-4)
ggsave(gg, filename = "__boxplot.beagle.tru.png", height = 10, width = 12-4)

ggsave(gg, filename = "__boxplot.beagle.err.pdf", height = 10, width = 12-4)
ggsave(gg, filename = "__boxplot.beagle.err.png", height = 10, width = 12-4)









##### 
stop()



by(d, list(d$type, d$fk), function(x) ks.exp.test(x$tru)$p.value)










x = by(p, list(p$type, p$fk), function(x) { rmsle(x$tru, x$det) })
x = t(array(x, dim(x), dimnames(x)))


ggplot(p) + facet_wrap(~type, ncol = 1) + stat_ecdf(aes(map, color = fk))

p = lapply(split(d, d$type), function(x) {
	z = density(x$map, bw = 0.01, n = 42, from=0, to=2.1)
	z = data.table(fk=x$fk[1], type=x$type[1], x = z$x, y = z$y / sum(z$y))
	z
})
p = rbindlist(p)

hist = ggplot(data=p) +
	facet_wrap(~type, ncol = 1) +
	#geom_point(aes(x=x, y=y), alpha=0.75, size=1) +
	geom_bar(aes(x=x, y=y), stat = "identity", colour = "white", size = 0.5, fill = "grey40") +
	#geom_histogram(aes(map), fill = "grey20") +
	geom_vline(xintercept=c(1), colour="white", size = 0.9) +
	geom_vline(xintercept=c(1), colour="black", size = 0.5, linetype = "22") +
	#geom_hline(yintercept = 0) +
	scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels =  c(0, 0.5, 1, 1.5, 2)) +
	#scale_fill_manual(values = c(DHG="limegreen", FGT="sienna2", HMM="royalblue2")) +
	coord_cartesian(ylim = c(0, 0.575), expand = F) +
	theme_few() +
	theme(aspect.ratio=1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = "top",
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(colour = "grey70"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	xlab("Relative distance to true breakpoint") +
	ylab("Density")
hist

ggsave(filename = sprintf("_plot.true-inferred.hist.%s.pdf", prefix), plot = hist, width = 9, height = 9)
### OR:
ggsave(filename = sprintf("_plot.beagle.hist.%s.pdf", prefix), plot = hist, width = 9, height = 6)



