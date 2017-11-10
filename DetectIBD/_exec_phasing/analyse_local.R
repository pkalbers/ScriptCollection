

infer.hap = function(tar, oth) {
	(tar + oth) / 4
}


files = dir(pattern = "^phase\\.rawdata\\.(.+)\\.RData$", path = "../rawdata")

for (file in sample(files)) {
	x = sub("^phase\\.rawdata\\.(.+)\\.RData$", "\\1", file)

	flagfile = sprintf("_flag.local.%s.txt", x)
	savefile = sprintf("phase.local.%s.RData", x)

	if (file.exists(savefile)) next

	Sys.sleep(runif(1, min = 0.01, max = 1))

	if (file.exists(flagfile)) next

	cat(".", file = flagfile, append = F)

	load(sprintf("../rawdata/%s", file))

	if (idv != as.numeric(x)) stop(sprintf("%d %s", idv, x))

	cat("\n\n", idv, "\n\n")


coord = cbind(this = idv, coord)


#
# local, estimated chunks
#

d = coord
d$underest = F
d$len = NA
d$het = NA
d$hht = NA
d$inf = NA
d$err = NA
d$err.foc = NA
d$err.sin = NA
d$err.frq = ""

for (i in 1:nrow(coord)) {
	if (i %% 50 == 0) cat(".")
	if (i %% 5000 == 0) cat("\n")

	if (coord$beg[i] < coord$beg[i] || coord$end[i] > coord$end[i]) {
		d$underest[i] = T
	}

	rng = (coord$beg[i]):(coord$end[i])

	tar = target[rng]
	oth = chunk[[ coord$key[i] ]]

	het = which(tar == 1)
	hht = intersect(het, which(oth == 1))
	inf = setdiff(het, hht)

	tru = NULL
	if (coord$h0[i]) tru = target.h0[rng]
	if (coord$h1[i]) tru = target.h1[rng]

	est = infer.hap(tar, oth)

	x = which(est == 0.5); if (length(x) > 0) est[x] = NA

	err = which(tru != round(est))

	foc = c()
	sin = c()

	err.frq = ""

	if (length(err) > 0) {
		foc = which(rng %in% unique(coord$foc))
		foc = which(err %in% foc)

		sin = which(rng %in% single)
		sin = which(err %in% sin)

		err.frq = paste(frq[rng][err], collapse = ",")
	}

	d$len[i] = length(rng)
	d$het[i] = length(het)
	d$hht[i] = length(hht)
	d$inf[i] = length(inf)
	d$err[i] = length(err)
	d$err.foc[i] = length(foc)
	d$err.sin[i] = length(sin)
	d$err.frq[i] = err.frq
}
cat("\n")

est.segment = d



#
# local, true chunks
#

d = coord
d$underest = F
d$len = NA
d$het = NA
d$hht = NA
d$inf = NA
d$err = NA
d$err.foc = NA
d$err.sin = NA
d$err.frq = ""

for (i in 1:nrow(coord)) {
	if (i %% 50 == 0) cat(".")
	if (i %% 5000 == 0) cat("\n")

	if (coord$true.beg[i] < coord$beg[i] || coord$true.end[i] > coord$end[i]) {
		d$underest[i] = T
		next
	}

	rng = (coord$true.beg[i]):(coord$true.end[i])

	tar = target[rng]
	oth = chunk[[ coord$key[i] ]]

	tmp = (coord$beg[i]):(coord$end[i])
	del = which(is.na(match(tmp, rng)))
	if (length(del) > 0) oth = oth[-del]

	het = which(tar == 1)
	hht = intersect(het, which(oth == 1))
	inf = setdiff(het, hht)

	tru = NULL
	if (coord$h0[i]) tru = target.h0[rng]
	if (coord$h1[i]) tru = target.h1[rng]

	est = infer.hap(tar, oth)

	x = which(est == 0.5); if (length(x) > 0) est[x] = NA

	err = which(tru != round(est))

	foc = c()
	sin = c()

	err.frq = ""

	if (length(err) > 0) {
		foc = which(rng %in% unique(coord$foc))
		foc = which(err %in% foc)

		sin = which(rng %in% single)
		sin = which(err %in% sin)

		err.frq = paste(frq[rng][err], collapse = ",")
	}

	d$len[i] = length(rng)
	d$het[i] = length(het)
	d$hht[i] = length(hht)
	d$inf[i] = length(inf)
	d$err[i] = length(err)
	d$err.foc[i] = length(foc)
	d$err.sin[i] = length(sin)
	d$err.frq[i] = err.frq
}
cat("\n")

tru.segment = d






#
# cluster, estimated chunks
#

cluster = split(coord, coord$foc)

d = NULL

k = 0
for (tag in names(cluster)) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")

	cls = cluster[[tag]]

	full = (min(cls$beg)):(max(cls$end))

	est0 = rep(0, length(full))
	est1 = rep(0, length(full))

	est = rep(NA, length(full))

	ncls = 0
	for (i in 1:nrow(cls)) {
		#if (cls$beg[i] < cls$beg[i] || cls$end[i] > cls$end[i]) next
		ncls = ncls + 1

		rng = (cls$beg[i]):(cls$end[i])
		rel = (rng - min(cls$beg) + 1)

		tar = target[rng]
		oth = chunk[[ cls$key[i] ]]

		tmp = infer.hap(tar, oth)

		x = which(tmp <  0.5); if (length(x) > 0) est0[rel[x]] = est0[rel[x]] + 1
		x = which(tmp >  0.5); if (length(x) > 0) est1[rel[x]] = est1[rel[x]] + 1
	}

	if (ncls == 0) next

	x = which(est0 > est1); if (length(x) > 0) est[x] = 0
	x = which(est1 > est0); if (length(x) > 0) est[x] = 1

	tar = target[full]
	het = which(tar == 1)
	hht = intersect(het, which(is.na(est)))
	inf = setdiff(het, hht)

	tru0 = target.h0[full]
	tru1 = target.h1[full]

	err = c()
	if (names(which.max(table(cls$h0)))[1] == "TRUE") {
		err = which(tru0 != est)
	} else {
		err = which(tru1 != est)
	}

	foc = c()
	sin = c()

	err.frq = ""

	if (length(err) > 0) {
		foc = which(full %in% unique(coord$foc))
		foc = which(err %in% foc)

		sin = which(full %in% single)
		sin = which(err %in% sin)

		err.frq = paste(frq[full][err], collapse = ",")
	}

	d = rbind(d, data.frame(fk = cls$fk[1],
													num = ncls,
													foc = cls$foc[1],
													len = length(full),
													het = length(het),
													hht = length(hht),
													inf = length(inf),
													err = length(err),
													err.foc = length(foc),
													err.sin = length(sin),
													err.frq = err.frq,
													stringsAsFactors = F))
}
cat("\n")

est.cluster = d


#
# cluster, true chunks
#

cluster = split(coord, coord$foc)

d = NULL

k = 0
for (tag in names(cluster)) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")

	cls = cluster[[tag]]

	full = (min(cls$true.beg)):(max(cls$true.end))

	est0 = rep(0, length(full))
	est1 = rep(0, length(full))

	est = rep(NA, length(full))

	ncls = 0
	for (i in 1:nrow(cls)) {
		if (cls$true.beg[i] < cls$beg[i] || cls$true.end[i] > cls$end[i]) next
		ncls = ncls + 1

		rng = (cls$true.beg[i]):(cls$true.end[i])
		rel = (rng - min(cls$true.beg) + 1)

		tar = target[rng]
		oth = chunk[[ cls$key[i] ]]

		tmp = (cls$beg[i]):(cls$end[i])
		del = which(is.na(match(tmp, rng)))
		if (length(del) > 0) oth = oth[-del]

		tmp = infer.hap(tar, oth)

		x = which(tmp <  0.5); if (length(x) > 0) est0[rel[x]] = est0[rel[x]] + 1
		x = which(tmp >  0.5); if (length(x) > 0) est1[rel[x]] = est1[rel[x]] + 1
	}

	if (ncls == 0) next

	x = which(est0 > est1); if (length(x) > 0) est[x] = 0
	x = which(est1 > est0); if (length(x) > 0) est[x] = 1

	tar = target[full]
	het = which(tar == 1)
	hht = intersect(het, which(is.na(est)))
	inf = setdiff(het, hht)

	tru0 = target.h0[full]
	tru1 = target.h1[full]

	err = c()
	if (names(which.max(table(cls$h0)))[1] == "TRUE") {
		err = which(tru0 != est)
	} else {
		err = which(tru1 != est)
	}

	foc = c()
	sin = c()

	err.frq = ""

	if (length(err) > 0) {
		foc = which(full %in% unique(coord$foc))
		foc = which(err %in% foc)

		sin = which(full %in% single)
		sin = which(err %in% sin)

		err.frq = paste(frq[full][err], collapse = ",")
	}

	d = rbind(d, data.frame(fk = cls$fk[1],
													num = ncls,
													foc = cls$foc[1],
													len = length(full),
													het = length(het),
													hht = length(hht),
													inf = length(inf),
													err = length(err),
													err.foc = length(foc),
													err.sin = length(sin),
													err.frq = err.frq,
													stringsAsFactors = F))
}
cat("\n")

tru.cluster = d


if (file.exists(savefile)) next

save(est.segment, tru.segment,
		 est.cluster, tru.cluster, file = savefile)

}




stop()

###

files = dir(pattern = "^phase\\.local\\.(.+)\\.RData$", path = "./local")

es = list()
ts = list()

ec = list()
tc = list()

k = 0
for (file in files) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")

	load(file)

	tag = sub("^phase\\.local\\.(.+)\\.RData$", "\\1", file)

	es[[tag]] = est.segment
	ts[[tag]] = tru.segment

	ec[[tag]] = est.cluster
	tc[[tag]] = tru.cluster
}


est.segment = es
tru.segment = ts

est.cluster = ec
tru.cluster = tc

save(est.segment, tru.segment,
		 est.cluster, tru.cluster, file = "phase.local.ALL.RData")



stop()

###

library(data.table)
library(ggplot2)
library(ggthemes)

load("phase.local.ALL.RData")

se <- function(x) {
	sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
}


###
### per individual
###

### frequencies of error sites

de = est.segment
dt = tru.segment

de = lapply(de, function(d) {
	k = sprintf("%d %d %d %d", d$this, d$other, d$beg, d$end)
	del = which(duplicated(k))
	if (length(del) > 0) {
		d = d[-del, ]
	}
	p = as.numeric(unlist(strsplit(d$err.frq, ",", T)))
	brk = seq(0, 5000, length.out = 201)
	del = which(p == 1); if (length(del) > 0) p = p[-del]
	p = cut(p, brk, include.lowest = T)
	p = table(p)
	data.table(x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
						 y = as.vector(p) / sum(p))
})

dt = lapply(dt, function(d) {
	k = sprintf("%d %d %d %d", d$this, d$other, d$true.beg, d$true.end)
	del = which(duplicated(k))
	if (length(del) > 0) {
		d = d[-del, ]
	}
	p = as.numeric(unlist(strsplit(d$err.frq, ",", T)))
	brk = seq(0, 5000, length.out = 201)
	del = which(p == 1); if (length(del) > 0) p = p[-del]
	p = cut(p, brk, include.lowest = T)
	p = table(p)
	data.table(x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
						 y = as.vector(p) / sum(p))
})


pe = data.table(type = "  Detected  ",
								x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
								mn = rep(0, 200),
								se = rep(0, 200),
								sd = rep(0, 200))
for (i in 1:200) {
	cat(".")
	tmp = sapply(de, function(x, i) x$y[i], i) * 100
	pe$mn[i] = mean(tmp)
	pe$se[i] = se(tmp)
	pe$sd[i] = sd(tmp)
}
cat("\n")

pt = data.table(type = "  True  ",
								x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
								mn = rep(0, 200),
								se = rep(0, 200),
								sd = rep(0, 200))
for (i in 1:200) {
	cat(".")
	tmp = sapply(dt, function(x, i) x$y[i], i) * 100
	pt$mn[i] = mean(tmp)
	pt$se[i] = se(tmp)
	pt$sd[i] = sd(tmp)
}
cat("\n")


p = rbind(pt, pe)
p$type = factor(p$type, levels = c("  True  ", "  Detected  "), ordered = T)
del = which(p$mn == 0)
if (length(del) > 0) {
	#p = p[-del, ]
	p$mn[del] = 1e-10
}


gg = ggplot(p) +
	geom_ribbon(aes(x=x, ymin=mn-se, ymax=mn+se, fill = type, colour = type), size = 0.5, alpha = 0.5) +
	geom_point(aes(x=x, y=mn, colour = type), size = 1.5) +
	scale_x_continuous(expand = c(0.02, 0)) +
	scale_y_log10(expand = c(0.02, 0), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = c(0.001, 0.01, 0.1, 1, 10, 100),
								minor_breaks = c(seq(0.001, 0.009, by=0.001), seq(0.01, 0.09, by=0.01), seq(0.1, 0.9, by=0.1), seq(1, 9, by=1), seq(10, 100, by=10))) +
	scale_fill_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	scale_colour_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	coord_cartesian(ylim = c(1e-3, 100)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		legend.justification = c(1, 1), legend.position = c(1, 1),
		legend.margin = unit(0.1, "cm"), legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		#legend.position = "bottom",
		legend.direction = "vertical",
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele frequency") + ylab("Mean density ± SE (%)")
gg
ggsave(gg, filename = "_plot.phase.errorsitefrqmnse.pdf", height = 4.5, width = 9)



###
### per fk
###

### segment error by fk

de = est.segment
dt = tru.segment

de = lapply(de, as.data.table)
de = rbindlist(de)
de = split(de, de$fk)

dt = lapply(dt, as.data.table)
dt = rbindlist(dt)
dt = split(dt, dt$fk)

# de = lapply(de, function(d) {
# 	k = sprintf("%d %d %d %d", d$this, d$other, d$beg, d$end)
# 	del = which(duplicated(k))
# 	if (length(del) > 0) {
# 		d = d[-del, ]
# 	}
# 	d
# })
#
# dt = lapply(dt, function(d) {
# 	k = sprintf("%d %d %d %d", d$this, d$other, d$true.beg, d$true.end)
# 	del = which(duplicated(k))
# 	if (length(del) > 0) {
# 		d = d[-del, ]
# 	}
# 	d
# })


getstats = function(d) {
	data.table(fk = d$fk[1],
						 het = mean(d$het / d$len, na.rm = T),
						 hht = mean(d$hht / d$len, na.rm = T),
						 inf = mean(d$inf / d$len, na.rm = T),
						 err = mean(d$err / d$len, na.rm = T),
						 foc = mean((d$err - d$err.foc) / d$len, na.rm = T),
						 sin = mean((d$err - d$err.sin) / d$len, na.rm = T),
						 focsin = mean((d$err - (d$err.foc + d$err.sin)) / d$len, na.rm = T),
						 se.het = se(d$het / d$len),
						 se.hht = se(d$hht / d$len),
						 se.inf = se(d$inf / d$len),
						 se.err = se(d$err / d$len),
						 se.foc = se((d$err - d$err.foc) / d$len),
						 se.sin = se((d$err - d$err.sin) / d$len),
						 se.focsin = se((d$err - (d$err.foc + d$err.sin)) / d$len))
}

pe = rbindlist(lapply(de, getstats))
pt = rbindlist(lapply(dt, getstats))

pe$type = "  Detected  "
pt$type = "  True  "

p = rbind(pt, pe)

p$type = factor(p$type, levels = c("  True  ", "  Detected  "), ordered = T)

### cluster error by fk

de = est.cluster
dt = tru.cluster

de = lapply(de, as.data.table)
de = rbindlist(de)
de = split(de, de$fk)

dt = lapply(dt, as.data.table)
dt = rbindlist(dt)
dt = split(dt, dt$fk)

pe = rbindlist(lapply(de, getstats))
pt = rbindlist(lapply(dt, getstats))

pe$type = "  Detected  "
pt$type = "  True  "

q = rbind(pt, pe)

q$type = factor(q$type, levels = c("  True  ", "  Detected  "), ordered = T)

p$mode = " Single segment"
q$mode = "Segment cluster"

d = rbind(p, q)

gg = ggplot(d) +
	facet_grid(.~mode) +
	geom_ribbon(aes(x=fk, ymin=sin-se.sin, ymax=sin+se.sin, fill = type, colour = type), size = 0.5, alpha = 0.5) +
	geom_point(aes(x=fk, y=sin, colour = type), size = 1.5, show.legend = F) +
	geom_ribbon(aes(x=fk, ymin=focsin-se.focsin, ymax=focsin+se.focsin, fill = type), size = 0.5, alpha = 0.5) +
	geom_point(aes(x=fk, y=focsin, colour = type), size = 1.5, shape = 1, show.legend = F) +
	scale_x_continuous(expand = c(0.02, 0), breaks = c(2, 5, 10, 15, 20, 25)) +
	scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 1, by = 0.02), minor_breaks = seq(0, 1, by = 0.01), labels = seq(0, 1, by = 0.02)*100) +
	# scale_y_log10(expand = c(0.02, 0), breaks = c(0.001, 0.01, 0.1, 1, 10), labels = c(0.001, 0.01, 0.1, 1, 10),
	# 							minor_breaks = c(seq(0.001, 0.009, by=0.001), seq(0.01, 0.09, by=0.01), seq(0.1, 0.9, by=0.1), seq(1, 9, by=1), seq(10, 100, by=10))) +
	scale_fill_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	scale_colour_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	#coord_cartesian(ylim = c(1e-3, 0.5)) +
	coord_cartesian(ylim = c(0, 0.15)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		legend.justification = c(1, 1), legend.position = c(1, 1)-0.01,
		legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.key.size = unit(0.5, "cm"),
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Mean error rate ± SE (%)")
gg
ggsave(gg, filename = "_plot.phase.fk-error.pdf", height = 4.5, width = 9)


###
### overall
###

### error site frqs

de = est.segment
dt = tru.segment

de = lapply(de, as.data.table)
de = rbindlist(de)

dt = lapply(dt, as.data.table)
dt = rbindlist(dt)

ke = sprintf("%d %d %d %d", de$this, de$other, de$beg, de$end)
kt = sprintf("%d %d %d %d", dt$this, dt$other, dt$true.beg, dt$true.end)

del = which(duplicated(ke))
if (length(del) > 0) {
	de = de[-del, ]
}
del = which(duplicated(kt))
if (length(del) > 0) {
	dt = dt[-del, ]
}


pe = as.numeric(unlist(strsplit(de$err.frq, ",", T)))
pt = as.numeric(unlist(strsplit(dt$err.frq, ",", T)))

brk = seq(0, 5000, length.out = 201)

del = which(pe == 1); if (length(del) > 0) pe = pe[-del]
del = which(pt == 1); if (length(del) > 0) pt = pt[-del]

pe = cut(pe, brk, include.lowest = T)
pe = table(pe)

pt = cut(pt, brk, include.lowest = T)
pt = table(pt)


p = rbind(data.table(type = "  True  ",
										 x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
										 y = (as.vector(pt) / sum(pt)) * 100),
					data.table(type = "  Detected  ",
										 x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
										 y = (as.vector(pe) / sum(pe)) * 100))

p$type = factor(p$type, levels = c("  True  ", "  Detected  "), ordered = T)

del = which(p$y == 0)
if (length(del) > 0) {
	#p = p[-del, ]
	p$y[del] = 1e-100
}



gg = ggplot(p) +
	geom_linerange(aes(x=x, ymax=y, ymin=1e-100, colour = type, fill = type), size = 1.4, alpha = 1/4, show.legend = F) +
	geom_point(aes(x=x, y=y, colour = type), shape = 15, size = 1.4) +
	scale_x_continuous(expand = c(0.003, 0), breaks = c(0, seq(0.1, 0.9, by=0.1), 1), labels = c("0", seq(0.1, 0.9, by=0.1), "1")) +
	scale_y_log10(expand = c(0.02, 0), breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), labels = c("0.0001", 0.001, 0.01, 0.1, 1, 10, 100),
								minor_breaks = c(seq(0.0001, 0.0009, by=0.0001), seq(0.001, 0.009, by=0.001), seq(0.01, 0.09, by=0.01), seq(0.1, 0.9, by=0.1), seq(1, 9, by=1), seq(10, 100, by=10))) +
	scale_fill_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	scale_colour_manual(values = c("  True  " = "black", "  Detected  " = "goldenrod")) +
	coord_cartesian(ylim = c(1e-3, 100)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		legend.justification = c(1, 1), legend.position = c(1, 1),
		legend.margin = unit(0.1, "cm"), legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		#legend.position = "bottom",
		legend.direction = "vertical",
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele frequency") + ylab("Density (%)")
gg
ggsave(gg, filename = "_plot.phase.errorsitefrq.pdf", height = 4.5, width = 9)


# stats

sfs = frq
del = which(sfs == 1); if (length(del) > 0) sfs = sfs[-del]
sfs = cut(sfs, brk, include.lowest = T)
sfs = table(sfs)
sfs = data.table(x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
								 y = (as.vector(sfs) / sum(sfs)) * 100)

spe = data.table(x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
								 y = (as.vector(pt) / sum(pt)) * 100)

spt = data.table(x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
								 y = (as.vector(pe) / sum(pe)) * 100)

cor(sfs$y, spe$y, method = "p")
cor(sfs$y, spe$y, method = "s")

cor(sfs$y, spt$y, method = "p")
cor(sfs$y, spt$y, method = "s")


### heterozygote distr

de = est.segment
dt = tru.segment

de = lapply(de, as.data.table)
de = rbindlist(de)

dt = lapply(dt, as.data.table)
dt = rbindlist(dt)

de = split(de, de$fk)
dt = split(dt, dt$fk)

gethets = function(d) {
	if (nrow(d) > 10000) { d = d[sample(1:nrow(d), 10000), ] }
	rbind(data.table(z = "Indeterminate   ", fk = d$fk[1],  mn = mean(d$hht / d$het, na.rm = T),           se = se(d$hht / d$het)),
				data.table(z = "Correct   ",      fk = d$fk[1],  mn = mean((d$inf - d$err) / d$het, na.rm = T), se = se((d$inf - d$err) / d$het)),
				data.table(z = "Incorrect   ",    fk = d$fk[1],  mn = mean(d$err / d$het, na.rm = T),           se = se(d$err / d$het)))
}

pe = rbindlist(lapply(de, gethets))
pt = rbindlist(lapply(dt, gethets))

pe$type = "  Detected  "
pt$type = "  True  "

p = rbind(pt, pe)

p$type = factor(p$type, levels = c("  True  ", "  Detected  "), ordered = T)


de = est.cluster
dt = tru.cluster

de = lapply(de, as.data.table)
de = rbindlist(de)

dt = lapply(dt, as.data.table)
dt = rbindlist(dt)

de = split(de, de$fk)
dt = split(dt, dt$fk)

pe = rbindlist(lapply(de, gethets))
pt = rbindlist(lapply(dt, gethets))

pe$type = "  Detected  "
pt$type = "  True  "

q = rbind(pt, pe)

q$type = factor(q$type, levels = c("  True  ", "  Detected  "), ordered = T)


p$mode = " Single segment"
q$mode = "Segment cluster"

d = rbind(p, q)

ggplot(d) +
	facet_wrap(type~mode, nrow = 2) +
	geom_ribbon(aes(x=fk, ymin=mn-se, ymax=mn+se, fill = z, colour = z), size = 0.5, alpha = 0.5) +
	geom_point(aes(x = fk, y = mn, colour = z), size = 0.75, show.legend = F) +
	scale_color_manual(values = c("Correct   " = "limegreen", "Incorrect   " = "red", "Indeterminate   " = "grey50")) +
	scale_fill_manual(values =  c("Correct   " = "limegreen", "Incorrect   " = "red", "Indeterminate   " = "grey50")) +
	scale_x_continuous(expand = c(0.02, 0), breaks = c(2, 5, 10, 15, 20, 25)) +
	scale_y_continuous(breaks = seq(0, 1, by=0.2), labels = seq(0, 1, by=0.2)*100) +
	coord_cartesian(ylim = c(0,1)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		legend.position = "bottom",
		legend.key.size = unit(0.5, "cm"),
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Mean proportion ± SE (%)")






stop()

data = est.segment
data = lapply(data, as.data.table)
data = rbindlist(data)

key = sprintf("%d %d %d %d", data$this, data$other, data$beg, data$end)
#key = sprintf("%d %d %d %d", data$this, data$other, data$true.beg, data$true.end)
del = which(duplicated(key))
if (length(del) > 0) {
	data = data[-del, ]
}


d = split(data, data$fk)

tmp = lapply(d, function(x) {
	as.numeric(unlist(strsplit(x$err.frq, ",", T)))
})
f = NULL
for (tag in names(tmp)[c(1,24)]) {
	z = tmp[[tag]]
	del = which(z == 1); if (length(del) > 0) z = z[-del]
	brk = seq(0, 5000, length.out = 201)
	z = cut(z, brk, include.lowest = T)
	z = table(z)
	f = rbind(f, data.table(fk = as.numeric(tag),
													brk = names(z),
													x = seq(25/5000/2, 1 - (25/5000/2), length.out = 200),
													y = as.vector(z) / sum(z)))
}

ggplot(f) + geom_point(aes(x=x, y=y, colour = factor(fk))) + scale_y_log10()



d = split(data, data$fk)


d = lapply(d, function(x) {
	data.table(fk = x$fk[1],
						 het = x$het / x$len,
						 hht = x$hht / x$het,
						 inf = x$inf / x$het,
						 err = x$err / x$het,
						 err.sin = x$err.sin / x$het,
						 err.nosin = (x$err - x$err.sin) / x$het)
})


a = lapply(d, function(d) {
	a = apply(d[, -(1)], 2, mean, na.rm = T)
	names(a) = names(d)[-(1)]
	a = as.data.table(t(a))
	cbind(fk = d$fk[1], a)
})
a = rbindlist(a)

b = lapply(d, function(d) {
	b = apply(d[, -(1)], 2, se)
	names(b) = names(d)[-(1)]
	b = as.data.table(t(b))
	cbind(fk = d$fk[1], b)
})
b = rbindlist(b)

names(b) = c("fk", sprintf("se.%s", names(b)[-1]))

p = cbind(a, b[, -1])


ggplot(p) +
	geom_ribbon(aes(x = fk, ymin = err - se.err, ymax = err + se.err), fill = "brown", alpha = 0.5) +
	geom_point(aes(x = fk, y = err), colour = "brown", alpha = 0.75) +
	geom_ribbon(aes(x = fk, ymin = err.nosin - se.err.nosin, ymax = err.nosin + se.err.nosin), fill = "purple1", alpha = 0.5) +
	geom_point(aes(x = fk, y = err.nosin), colour = "purple1", alpha = 0.75) +
	#scale_y_log10() +
	scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, by = 0.01)) +
	scale_x_continuous(breaks = c(2, 5, 10, 15, 20, 25)) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.border = element_rect(fill = NA, color = "black", size = 0.5),
				#panel.grid.y = element_line(colour = "grey80", size = 0.5),
				panel.grid.major.y = element_line(colour = "grey80", size = 0.25),
				panel.grid.minor.y = element_blank())









apply(d, 2, mean, na.rm = T)
apply(d, 2, se)



d = lapply(tru.segment, function(x) {
	x = lapply(split(x, x$fk), function(x) {
		data.table(this = x$this[1],
							 fk = x$fk[1],
							 het = x$het / x$len,
							 hht = x$hht / x$het,
							 inf = x$inf / x$het,
							 err = x$err / x$het,
							 err.sin = x$err.sin / x$het)
	})
	rbindlist(x)
})
d = rbindlist(d)
d = split(d, d$fk)

mn = lapply(d, function(d) {
	mn = apply(d[, -(1:2)], 2, mean, na.rm = T)
	names(mn) = names(d)[-(1:2)]
	mn = as.data.table(t(mn))
	cbind(fk = d$fk[1], mn)
})
mn = rbindlist(mn)

se = lapply(d, function(d) {
	se = apply(d[, -(1:2)], 2, se)
	names(se) = names(d)[-(1:2)]
	se = as.data.table(t(se))
	cbind(fk = d$fk[1], se)
})
se = rbindlist(se)

names(se) = c("fk", sprintf("se.%s", names(se)[-1]))

p = cbind(mn, se[, -1])

ggplot(p) +
	geom_point(aes(x = fk, y = err)) +
	geom_line(aes(x = fk, y = err - se.err)) +
	geom_line(aes(x = fk, y = err + se.err)) +
	scale_y_log10()












