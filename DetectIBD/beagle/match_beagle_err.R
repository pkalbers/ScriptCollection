

library(data.table)


load("../result.truth.local.RData")


load("data.OutOfAfricaHapMap20.GenErr_1000G.RData")

pos = round(POS)

i = anyDuplicated(pos)
while (i != 0) {
	if (i == 1) {
		pos[i] = pos[i] - 1
	} else {
		pos[i] = pos[i] + 1
	}
	i = anyDuplicated(pos)
}

truth$position = pos[match(truth$position, POS)]
truth$lhs = pos[truth$lhs.index]
truth$rhs = pos[truth$rhs.index]


tru = as.data.table(truth)



full = tru[order(tru$fk),]


key = sprintf("%d %d %d %d", tru$h0, tru$h1, tru$lhs.index, tru$rhs.index)
del = which(duplicated(key))
tru = tru[-del,]

uniq = tru[order(tru$fk),]



full$pair = sprintf("%d %d", full$h0, full$h1)
uniq$pair = sprintf("%d %d", uniq$h0, uniq$h1)





d = fread("./beagle_ibd.OutOfAfricaHapMap20.GenErr_1000G.H.ibd", header = F, stringsAsFactors = F)
d = fread("./beagle_ibd.OutOfAfricaHapMap20.GenErr_1000G.P.ibd", header = F, stringsAsFactors = F)


names(d) = c("indv1", "x1", "indv2", "x2", "chr", "beg", "end", "lod")

d$g0 = as.numeric(sub("^ID([0-9]+)$", "\\1", d$indv1))
d$g1 = as.numeric(sub("^ID([0-9]+)$", "\\1", d$indv2))

d$h0 = (d$g0 * 2) + (d$x1 - 2)
d$h1 = (d$g1 * 2) + (d$x2 - 2)

d$tag = 1:nrow(d)

d$pair = sprintf("%d %d", d$h0, d$h1)


dh = d
dp = d

# dp1 = dp
# dp2 = dp
# dp3 = dp
#
# z = (dp$h0 %% 2 == 0)
# x = which(z)
# y = which(!z)
# dp1$h0[x] = dp1$h0[x] - 1
# dp1$h0[y] = dp1$h0[y] + 1
# dp2$h0[x] = dp2$h0[x] - 1
# dp2$h0[y] = dp2$h0[y] + 1
#
# z = (dp$h1 %% 2 == 0)
# x = which(z)
# y = which(!z)
# dp2$h1[x] = dp2$h1[x] - 1
# dp2$h1[y] = dp2$h1[y] + 1
# dp3$h1[x] = dp3$h1[x] - 1
# dp3$h1[y] = dp3$h1[y] + 1
#
# dp1$pair = sprintf("%d %d", dp1$h0, dp1$h1)
# dp2$pair = sprintf("%d %d", dp2$h0, dp2$h1)
# dp3$pair = sprintf("%d %d", dp3$h0, dp3$h1)
#
# dg = dp
# dg$pair = sprintf("%d %d", dg$g0, dg$g1)


length(unique(dh$pair))  #  5807497
length(unique(dp$pair))  #  5937559

# length(unique(dp1$pair))  #
# length(unique(dp2$pair))  #
# length(unique(dp3$pair))  #
# length(unique(dg$pair))  #

length(intersect(unique(dh$pair), unique(dp$pair)))  #  3258158

# length(intersect(unique(dp$pair), unique(dp1$pair)))  #
# length(intersect(unique(dp$pair), unique(dp2$pair)))  #
# length(intersect(unique(dp$pair), unique(dp3$pair)))  #

choose(5000, 2)  #  12497500
length(unique(full$pair))  #  2638412

median((dh$end - dh$beg) + 1)
median((dp$end - dp$beg) + 1)
median((uniq$rhs - uniq$lhs) + 1)


# # phasing correction
#
# load("data.OutOfAfricaHapMap20.H.RData")
# load("data.OutOfAfricaHapMap20.P.RData")
#
# hs = apply(H, 2, function(h) sum(h+1))
# ps = apply(P, 2, function(h) sum(h+1))
#
# swap = 1:length(hs)
# k = 0
# for (j in seq(2, length(hs), by = 2)) {
# 	i = j -1
#
# 	if (abs(hs[i] - ps[j]) < abs(hs[i] - ps[i]) & (abs(hs[j] - ps[i]) < abs(hs[j] - ps[j]))) {
# 		x = swap[i]
# 		swap[i] = swap[j]
# 		swap[j] = x
# 		k = k + 1
# 	}
# }




# segment overlap

length(unique(full$pair))  #  2638412
pairs.h = intersect(unique(dh$pair), unique(uniq$pair))  #  2332363
pairs.p = intersect(unique(dp$pair), unique(uniq$pair))  #  1661487

pairs = intersect(pairs.h, pairs.p)  #  1522900
pairs = sample(pairs, size = 10000)


d = dh
d = dp


dd = d[which(d$pair %in% pairs), ]
tt = uniq[which(uniq$pair %in% pairs), ]

t.dd = table(dd$pair)
t.tt = table(tt$pair)


ddd = dd[which(dd$pair %in% pairs),]
ttt = tt[which(tt$pair %in% pairs),]

ddd = ddd[order(ddd$pair),]
ttt = ttt[order(ttt$pair),]


beg = split(ddd$beg, ddd$pair)
end = split(ddd$end, ddd$pair)

lhs = split(ttt$lhs, ttt$pair)
rhs = split(ttt$rhs, ttt$pair)


p.tru = c()
p.est = c()

i = 0
for (pair in pairs) {
	i = i + 1
	if (i %% 1000 == 0) cat(i, "\n")

	xl = beg[[pair]]
	xr = end[[pair]]

	yl = lhs[[pair]]
	yr = rhs[[pair]]

	for (j in 1:length(xl)) {
		for (k in 1:length(yl)) {
			if (xl[j] <= yl[k]) {
				if (xr[j] >= yl[k]) {
					l = (min(c(yr[k], xr[j])) - max(c(yl[k], xl[j])) + 1)
					p.tru = c(p.tru, l / (yr[k] - yl[k] + 1))
					p.est = c(p.est, l / (xr[j] - xl[j] + 1))
				}
			} else if (xl[j] <= yr[k]) {
				l = (min(c(yr[k], xr[j])) - max(c(yl[k], xl[j])) + 1)
				p.tru = c(p.tru, l / (yr[k] - yl[k] + 1))
				p.est = c(p.est, l / (xr[j] - xl[j] + 1))
			}
		}
	}
}



mean(p.tru)  #  H : 0.1821816  P : 0.1828235
mean(p.est)  #  H : 0.9882891  P : 0.9865589


library(ggplot2)
library(ggthemes)


qh = rbind(data.table(p = p.tru, tag = "Relative to true segment",     mode = "(a) True haplotypes"),
					 data.table(p = p.est, tag = "Relative to inferred segment", mode = "(a) True haplotypes"))

qp = rbind(data.table(p = p.tru, tag = "Relative to true segment",     mode = "(b) Phased haplotypes"),
					 data.table(p = p.est, tag = "Relative to inferred segment", mode = "(b) Phased haplotypes"))

q = rbind(qh, qp)


ggplot(q) +
	facet_grid(.~mode) +
	stat_bin(aes(p*100, ..count.., color = tag), bins = 21, position = "identity", geom = "line", size=2/3) +
	stat_bin(aes(p*100, ..count.., color = tag, shape = tag), bins = 21, position = "identity", geom = "point", size=1.5) +
	scale_x_continuous(breaks = seq(0, 100, length.out = 11)) +
	scale_y_continuous(minor_breaks = seq(0,1e5, by = 1000)) +
	scale_color_ptol() +
	scale_shape_circlefill() +
	coord_cartesian(xlim = c(0,100)) +
	theme_few() +
	theme(aspect.ratio = 0.8, #16/9,
				panel.border = element_rect(fill = NA, colour = "black", size = 1/3),
				panel.background = element_rect(fill = "grey90", colour = "black", size = 1/3),
				panel.grid.minor = element_line(colour = "grey98", size = 1/3),
				panel.grid.major = element_line(colour = "grey98", size = 1/2),
				strip.text.x = element_text(face = "bold", hjust = 0),
				legend.position = c(1-0.99, 0.98), legend.justification = c(0, 1),
				legend.background = element_rect(fill = "white", colour = "grey50", size = 1/3),
				legend.margin = margin(0, 2, 1, 1, unit = "mm"),
				legend.title = element_blank()) +
	xlab("Overlap (%)") + ylab("Count")

ggsave(filename = "__plot.beagle_overlaps.pdf", width = 8, height = 4)








# exhaustive matching

d = dh
d = dp


pairs = intersect(unique(d$pair), unique(full$pair))

length(pairs)  #  H : 2089489 , P : 1661487

ddd = d[which(d$pair %in% pairs), ]
ttt = full[which(full$pair %in% pairs), ]


ddd = ddd[order(ddd$pair),]
ttt = ttt[order(ttt$pair),]


ddd$tag = 1:nrow(ddd)  #  H :  6589784  P : 3740239
ttt$tag = 1:nrow(ttt)  #  H : 10364016  P : 7538333


beg = split(ddd$beg, ddd$pair)
end = split(ddd$end, ddd$pair)
pos = split(ttt$position, ttt$pair)

xxx = split(ddd$tag, ddd$pair)
yyy = split(ttt$tag, ttt$pair)

zbeg = split(beg, cut(1:length(beg), breaks = 1000))
zend = split(end, cut(1:length(end), breaks = 1000))
zpos = split(pos, cut(1:length(pos), breaks = 1000))

zxxx = split(xxx, cut(1:length(xxx), breaks = 1000))
zyyy = split(yyy, cut(1:length(yyy), breaks = 1000))




qwerty = function(pair, beg, end, pos, xxx, yyy) {
	ll = beg[[pair]]
	rr = end[[pair]]
	pp = pos[[pair]]

	xx = xxx[[pair]]
	yy = yyy[[pair]]

	idx = c()

	for (k in 1:length(ll)) {

		hit = which(ll[k] <= pp & pp <= rr[k])

		if (length(hit) > 0) {

			idx = c(idx, sprintf("%d %d", xx[k], yy[hit]))

		}

	}
	idx
}


idx = list()

i = 0
for (tag in intersect(intersect(names(zbeg), names(zend)), names(zpos))) {
	i = i + 1
	if (i %% 10 == 0) cat(sprintf("%d (%d)\n", i, length(idx)))

	beg = zbeg[[tag]]
	end = zend[[tag]]
	pos = zpos[[tag]]

	xxx = zxxx[[tag]]
	yyy = zyyy[[tag]]


	idx = c(idx, unlist(lapply(intersect(intersect(names(beg), names(end)), names(pos)), qwerty, beg, end, pos, xxx, yyy)))
}


idx = unlist(idx)

length(idx)  #  H : 4590639  P : 1165436


xd = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", idx))
xt = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", idx))


range(ttt$h0[xt] - ddd$h0[xd])
range(ttt$h1[xt] - ddd$h1[xd])


all = data.table(position = ttt$position[xt],
								 fk = ttt$fk[xt],
								 h0 = ttt$h0[xt],
								 h1 = ttt$h1[xt],
								 wall = ttt$wall[xt],
								 t.lhs = ttt$lhs[xt],
								 t.rhs = ttt$rhs[xt],
								 b.lhs = ddd$beg[xd],
								 b.rhs = ddd$end[xd])

all = all[order(all$fk, sample(1:nrow(all))),]

key = sprintf("%d %d %d %d", all$h0, all$h1, all$t.lhs, all$t.rhs)
del = which(duplicated(key))

sub = all[-del,]

key = sprintf("%d %d %d %d", sub$h0, sub$h1, sub$b.lhs, sub$b.rhs)
del = which(duplicated(key))
if (length(del) > 0) {
	sub = sub[-del,]
}

nrow(sub)  #  H : 1504613  P : 395785



allH = all
subH = sub

allP = all
subP = sub



save(allH, allP, subH, subP, file = "match.beagle.all_sub_segments.RData")












