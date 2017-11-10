

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
