

library(data.table)



marker = fread("./truH.marker.txt", header = T)
times = fread("./OutOfAfricaHapMap20.times.txt", header = T)

d = cbind(marker, times[marker$MarkerID + 1, ])


k = which(d$AlleleCount1 >= 0 & d$AlleleCount1 <= 55)
d = d[k,]

k = which(d$AlleleCount1 == d$GenotypeCount1)
d = d[k,]


d$fk = d$AlleleCount1

tmp = fread("./errH.marker.txt", header = T)

d$ffk = tmp$AlleleCount1[ d$MarkerID + 1 ]


# d = split(d, d$fk)
# d = lapply(d, function(x) {
# 	if (nrow(x) < 5000) return(x)
# 	x[sample(1:nrow(x), 5000), ]
# })
# d = rbindlist(d)

table(d$fk)

d = d[order(d$MarkerID), ]


d$fTru = cut(d$fk, breaks = c(2,5,10,15,20,25,30,35,40,45,50,55), include.lowest = T)
d$fErr = cut(d$ffk, breaks = c(2,5,10,15,20,25,30,35,40,45,50,55), include.lowest = T)


load("~/Research/DetectIBD/data.OutOfAfricaHapMap20.H.RData")
tru = H
load("~/Research/DetectIBD/GenErr_1000G/data.OutOfAfricaHapMap20.GenErr_1000G.H.RData")
err = H
rm(H)
gc()


tru = tru[d$MarkerID + 1, ]
err = err[d$MarkerID + 1, ]

d$cTru = rowSums(tru)
d$cErr = rowSums(err)


xxx = tru - err
rm(tru)
rm(err)
gc()


d$fp = apply(xxx, 1, function(x) length(which(x < 0)) )
d$fn = apply(xxx, 1, function(x) length(which(x > 0)) )

rm(xxx)
gc()

del = which(d$fk != d$cTru)
del = which(d$ffk != d$cErr)


table(d$fk[ which(d$fp == 0 & d$fn == 0) ])
table(d$fk[ which(d$fp != 0 & d$fn == 0) ])
table(d$fk[ which(d$fp == 0 & d$fn != 0) ])


id = d[ which(d$fp == 0 & d$fn == 0), ]
fp = d[ which(d$fp != 0 & d$fn == 0), ]
fn = d[ which(d$fp == 0 & d$fn != 0), ]
ff = d[ which(d$fp != 0 & d$fn != 0), ]

table(id$fErr)
table(fp$fErr)
table(fn$fErr)
table(ff$fErr)


tmp.f = function(x) {
	if (nrow(x) < 1000) return(x)
	x[sample(1:nrow(x), 1000, replace = F, prob = x$cErr / (max(x$cErr) + 1)),]
}

id = rbindlist(lapply(split(id, id$fErr), tmp.f))
fp = rbindlist(lapply(split(fp, fp$fErr), tmp.f))
fn = rbindlist(lapply(split(fn, fn$fErr), tmp.f))
ff = rbindlist(lapply(split(ff, ff$fErr), tmp.f))

del = which(id$fErr == "(50,55]"); id = id[-del,]
del = which(fp$fErr == "(50,55]"); fp = fp[-del,]
del = which(fn$fErr == "(50,55]"); fn = fn[-del,]
del = which(ff$fErr == "(50,55]"); ff = ff[-del,]


x<-order(id$MarkerID); id = id[x,]
x<-order(fp$MarkerID); fp = fp[x,]
x<-order(fn$MarkerID); fn = fn[x,]
x<-order(ff$MarkerID); ff = ff[x,]


table(id$fk)
table(fp$fk)
table(fn$fk)
table(ff$fk)

table(id$ffk)
table(fp$ffk)
table(fn$ffk)
table(ff$ffk)


cat(id$Position, sep = "\n", file = "expos_id.txt")
cat(fn$Position, sep = "\n", file = "expos_fn.txt")
cat(fp$Position, sep = "\n", file = "expos_fp.txt")
cat(ff$Position, sep = "\n", file = "expos_ff.txt")

save(d, id, fp, fn, ff, file = "expos.id_fp_fn_ff.RData")
load("expos.id_fp_fn_ff.RData")



library(ggplot2)
library(ggthemes)

ggplot(d[which(d$fk <= 50 & d$fk == d$ffk),]) +
	geom_bar(aes(x = factor(fk), fill = factor(fp)), position = "fill", width = 1) +
	geom_vline(xintercept = 2:56 - 0.5, size = 1/3, color = "grey60") +
	scale_fill_brewer(type = "seq", palette = "YlOrRd") +
	scale_y_continuous(breaks = (1:9)/10) +
	coord_cartesian(expand = F) +
	theme_few() +
	theme(panel.border = element_rect(fill = NA, colour = "grey40", size = 1/2),
				legend.title = element_blank()) +
	xlab("Allele count") + ylab("Proportion") + ggtitle("False positives")
ggsave(filename = "_plot.ErrorDistr.FP_eq.pdf", width = 16*(3/4), height = 9*(3/4))

ggplot(d[which(d$fk <= 50 & d$fk == d$ffk),]) +
	geom_bar(aes(x = factor(fk), fill = factor(fn)), position = "fill", width = 1) +
	geom_vline(xintercept = 2:56 - 0.5, size = 1/3, color = "grey60") +
	scale_fill_brewer(type = "seq", palette = "YlOrRd") +
	scale_y_continuous(breaks = (1:9)/10) +
	coord_cartesian(expand = F) +
	theme_few() +
	theme(panel.border = element_rect(fill = NA, colour = "grey40", size = 1/2),
				legend.title = element_blank()) +
	xlab("Allele count") + ylab("Proportion") + ggtitle("False negatives")
ggsave(filename = "_plot.ErrorDistr.FN_eq.pdf", width = 16*(3/4), height = 9*(3/4))





stop()


d$age = exp(log(d$node.time) + ((log(d$parent.time) - log(d$node.time)) / 2))

d$int = (cut((d$age), breaks = 10^(seq(log10(1), log10(1e4), length.out = 5)), include.lowest = T))


expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}

d$exp = sapply(d$AlleleCount1, function(x) expected.age(x, 5000))


x = table(d$AlleleCount1, d$int)
print(x)
table(d$AlleleCount1)/nrow(d)
table(d$int)/nrow(d)

z = t(t(x) / rowSums(t(x)))


d$pr = NA

for (r in rownames(z)) {
	cat(r, ":")
	a = as.numeric(r)
	a = which(d$AlleleCount1 == a)
	for (c in colnames(z)) {
		cat(c, "")
		b = which(d$int == c)
		i = intersect(a, b)
		if (length(i) == 0) next
		d$pr[i] = z[r,c]
	}
	cat("\n")
}

del = which(is.na(d$pr))
if (length(del) > 0) d = d[-del,]

range(d$pr)

s = sample(1:nrow(d), size = 10000, replace = F, prob = d$pr^2)

table(d$AlleleCount1[s])
table(d$int[s])
table(d$AlleleCount1[s], d$int[s])


a = split(d, d$AlleleCount1)
for (tag in names(a)) {
	print(tag)
	a[[tag]]$a = NA
	sub = (1 - table(a[[tag]]$int) / nrow(a[[tag]]))
	for (x in names(sub)) {
		if (sub[x] == 1) next
		a[[tag]]$a[which(a[[tag]]$int == x)] = sub[x]
	}
}
a = rbindlist(a)
del = which(is.na(a$int))
if (length(del) > 0) a = a[-del,]

table(a$AlleleCount1)
table(a$int)


a = split(d, d$AlleleCount1)
a = lapply(a, function(x, mn) {
	bb = table(x$int)
	b = (1-(bb / sum(bb)))
	x$lpr = 0
	for (tag in names(b)) {
		if (bb[tag] == 0) next
		x$lpr[which(x$int == tag)] = b[tag]
	}
	x$pr = x$pr / sum(x$pr)
	x$lpr = x$lpr / sum(x$lpr)
	x[sample(1:nrow(x), mn, replace = F, prob = x$pr^2 * x$lpr^2),]
	#x[sample(1:nrow(x), mn, replace = F),]
}, 500 ) # min(table(d$AlleleCount1)))
a = rbindlist(a)
table(a$AlleleCount1)
table(a$int)
table(a$AlleleCount1, a$int)




b = split(d, d$int)
b = lapply(b, function(x, mn) {
	if (is.null(x)) return(NULL)
	x[sample(1:nrow(x), mn),]
}, min(table(d$int)))
b = rbindlist(b)
table(b$AlleleCount1)
table(b$int)



mt = fread("~/Research/rvage/africa/truH.marker.txt")
me = fread("~/Research/rvage/africa/errH.marker.txt")


i = which(mt$AlleleCount1 == mt$GenotypeCount1)

k0 = intersect(which(mt$AlleleCount1 > 1), which(me$AlleleCount1 > 1))
k1 = intersect(which(mt$AlleleCount1 <= 61), which(me$AlleleCount1 <= 61))
k = intersect(k0, k1)

k = intersect(k, i)


f.id = intersect(which(mt$AlleleCount1 == me$AlleleCount1), k)
f.fn = intersect(which(mt$AlleleCount1 >  me$AlleleCount1), k)
f.fp = intersect(which(mt$AlleleCount1 <  me$AlleleCount1), k)


tmp = split(f.id, me$AlleleCount1[f.id])
tmp = split(f.fn, me$AlleleCount1[f.fn])
tmp = split(f.fp, me$AlleleCount1[f.fp])


tmp = lapply(tmp, function(x) {
	if (length(x) < 100) return(x)
	sort(sample(x, 100))
})
tmp = unlist(tmp)

del = which(me$AlleleCount1[tmp] > 51)
if (length(del) > 0) tmp = tmp[-del]

table(me$AlleleCount1[tmp])


p.id = me$Position[tmp]
p.fn = me$Position[tmp]
p.fp = me$Position[tmp]


cat(p.id, sep = "\n", file = "~/Research/rvage/africa/expos.f_id.txt")
cat(p.fn, sep = "\n", file = "~/Research/rvage/africa/expos.f_fn.txt")
cat(p.fp, sep = "\n", file = "~/Research/rvage/africa/expos.f_fp.txt")
