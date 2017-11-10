

infer.hap = function(tar, oth) {
	(tar + oth) / 4
}


files = dir(pattern = "^phase\\.rawdata\\.(.+)\\.RData$", path = "../rawdata")

for (file in sample(files)) {
	x = sub("^phase\\.rawdata\\.(.+)\\.RData$", "\\1", file)

	flagfile = sprintf("_flag.shap.%s.txt", x)
	savefile = sprintf("phase.shap.%s.RData", x)

	if (file.exists(savefile)) next

	Sys.sleep(runif(1, min = 0.01, max = 1))

	if (file.exists(flagfile)) next

	cat(".", file = flagfile, append = F)

	load(sprintf("../rawdata/%s", file))

	if (idv != as.numeric(x)) stop(sprintf("%d %s", idv, x))

	cat("\n\n", idv, "\n\n")


	coord = cbind(this = idv, coord)


	#
	# cluster, estimated chunks
	#

	cluster = split(coord, coord$foc)

	shap = list()

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
		tru0 = target.h0[full]
		tru1 = target.h1[full]

		shap[[ tag ]] = list(fk = cls$fk[1],
												 is0 = cls$h0[1],
												 is1 = cls$h1[1],
												 this = idv,
												 others = cls$other,
												 rng = full,
												 tar = tar,
												 tru0 = tru0,
												 tru1 = tru1,
												 est = est)

	}
	cat("\n")

	est.shap = shap


	#
	# cluster, true chunks
	#

	cluster = split(coord, coord$foc)

	shap = list()

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
		tru0 = target.h0[full]
		tru1 = target.h1[full]

		shap[[ tag ]] = list(fk = cls$fk[1],
												 is0 = cls$h0[1],
												 is1 = cls$h1[1],
												 this = idv,
												 others = cls$other,
												 rng = full,
												 tar = tar,
												 tru0 = tru0,
												 tru1 = tru1,
												 est = est)

	}
	cat("\n")

	tru.shap = shap


	if (file.exists(savefile)) next

	save(est.shap, tru.shap, file = savefile)

}





stop()
#######


load("../frq.RData")


flip.hap = function(h) {
	x = which(!is.na(h))
	h[x] = abs(h[x] - 1)
	h
}


#
# estimated
#

sites = as.numeric(names(est.shap))

d = NULL

k = 0
for (tag in names(est.shap)) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")

	site = as.numeric(tag)
	shap = est.shap[[tag]]

	if (shap$is0 == shap$is1) next

	match = which(sites %in% shap$rng)
	if (length(match) == 0) next
	match = as.character(sites[match])

	est = rep(NA, length(shap$rng))

	est0 = rep(0, length(shap$rng))
	est1 = rep(0, length(shap$rng))

	x = which(shap$est == 0); if (length(x) > 0) { est0[x] = est0[x] + 1 }
	x = which(shap$est == 1); if (length(x) > 0) { est1[x] = est1[x] + 1 }

	n = 0
	nr = 0
	for (mat in match) {
		sit = as.numeric(mat)
		if (sit == site) next
		sub = est.shap[[mat]]

		if (sub$is0 == sub$is1) next

		ovr = intersect(shap$rng, sub$rng)
		rr = match(ovr, shap$rng)
		aa = match(ovr, sub$rng)

		tar = shap$tar[rr]
		het = which(tar == 1)

		ref = shap$est[rr]
		alt0 = sub$est[aa]
		alt1 = flip.hap(sub$est[aa])

		x = which(!is.na(ref[het]) & !is.na(alt0[het]))

		if (length(x) == 0) next

		n0 = length(which(ref[het][x] == alt0[het][x]))
		n1 = length(which(ref[het][x] == alt1[het][x]))

		if (n0 == n1) next
		if (n0 > n1) {
			x = which(alt0[het] == 0); if (length(x) > 0) { est0[rr][het][x] = est0[rr][het][x] + 1 }
			x = which(alt0[het] == 1); if (length(x) > 0) { est1[rr][het][x] = est1[rr][het][x] + 1 }
			if (shap$is0 == sub$is0) nr = nr + 1
		}
		if (n1 > n0) {
			x = which(alt1[het] == 0); if (length(x) > 0) { est0[rr][het][x] = est0[rr][het][x] + 1 }
			x = which(alt1[het] == 1); if (length(x) > 0) { est1[rr][het][x] = est1[rr][het][x] + 1 }
			if (shap$is0 != sub$is0) nr = nr + 1
		}

		n = n + 1
	}

	if (n == 0) next

	x = which(est0 > est1); if (length(x) > 0) est[x] = 0
	x = which(est1 > est0); if (length(x) > 0) est[x] = 1

	het = which(shap$tar == 1)
	hht = which(is.na(est))
	inf = intersect(which(!is.na(est)), het)

	tru = NULL
	if (shap$is0 && !shap$is1) {
		tru = shap$tru0
	} else if (!shap$is0 && shap$is1) {
		tru = shap$tru1
	} else {
		next
	}

	err = which(tru[inf] != est[inf])
	err.foc = c()
	err.sin = c()
	err.frq = c()

	if (length(err) > 0) {
		err.foc = intersect(sites, shap$rng[inf][err])
		err.frq = frq[shap$rng[inf][err]]
		err.sin = which(err.frq == 1)
	}

	tmp = data.frame(fk = shap$fk,
									 this = shap$this,
									 foc = site,
									 len = length(shap$rng),
									 num = n,
									 dis = nr,
									 het = length(het),
									 hht = length(hht),
									 inf = length(inf),
									 err = length(err),
									 err.foc = length(err.foc),
									 err.sin = length(err.sin),
									 err.frq = paste(err.frq, collapse = ","))

	d = rbind(d, tmp)
}

est.stack = d


#
# true
#

sites = as.numeric(names(tru.shap))

d = NULL

k = 0
for (tag in names(tru.shap)) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")

	site = as.numeric(tag)
	shap = tru.shap[[tag]]

	if (shap$is0 == shap$is1) next

	match = which(sites %in% shap$rng)
	if (length(match) == 0) next
	match = as.character(sites[match])

	est = rep(NA, length(shap$rng))

	est0 = rep(0, length(shap$rng))
	est1 = rep(0, length(shap$rng))

	x = which(shap$est == 0); if (length(x) > 0) { est0[x] = est0[x] + 1 }
	x = which(shap$est == 1); if (length(x) > 0) { est1[x] = est1[x] + 1 }

	n = 0
	nr = 0
	for (mat in match) {
		sit = as.numeric(mat)
		if (sit == site) next
		sub = tru.shap[[mat]]

		if (sub$is0 == sub$is1) next

		ovr = intersect(shap$rng, sub$rng)
		rr = match(ovr, shap$rng)
		aa = match(ovr, sub$rng)

		tar = shap$tar[rr]
		het = which(tar == 1)

		ref = shap$est[rr]
		alt0 = sub$est[aa]
		alt1 = flip.hap(sub$est[aa])

		x = which(!is.na(ref[het]) & !is.na(alt0[het]))

		if (length(x) == 0) next

		n0 = length(which(ref[het][x] == alt0[het][x]))
		n1 = length(which(ref[het][x] == alt1[het][x]))

		if (n0 == n1) next
		if (n0 > n1) {
			x = which(alt0[het] == 0); if (length(x) > 0) { est0[rr][het][x] = est0[rr][het][x] + 1 }
			x = which(alt0[het] == 1); if (length(x) > 0) { est1[rr][het][x] = est1[rr][het][x] + 1 }
			if (shap$is0 == sub$is0) nr = nr + 1
		}
		if (n1 > n0) {
			x = which(alt1[het] == 0); if (length(x) > 0) { est0[rr][het][x] = est0[rr][het][x] + 1 }
			x = which(alt1[het] == 1); if (length(x) > 0) { est1[rr][het][x] = est1[rr][het][x] + 1 }
			if (shap$is0 != sub$is0) nr = nr + 1
		}

		n = n + 1
	}

	if (n == 0) next

	x = which(est0 > est1); if (length(x) > 0) est[x] = 0
	x = which(est1 > est0); if (length(x) > 0) est[x] = 1

	het = which(shap$tar == 1)
	hht = which(is.na(est))
	inf = intersect(which(!is.na(est)), het)

	tru = NULL
	if (shap$is0 && !shap$is1) {
		tru = shap$tru0
	} else if (!shap$is0 && shap$is1) {
		tru = shap$tru1
	} else {
		next
	}

	err = which(tru[inf] != est[inf])
	err.foc = c()
	err.sin = c()
	err.frq = c()

	if (length(err) > 0) {
		err.foc = intersect(sites, shap$rng[inf][err])
		err.frq = frq[shap$rng[inf][err]]
		err.sin = which(err.frq == 1)
	}

	tmp = data.frame(fk = shap$fk,
									 this = shap$this,
									 foc = site,
									 len = length(shap$rng),
									 num = n,
									 dis = nr,
									 het = length(het),
									 hht = length(hht),
									 inf = length(inf),
									 err = length(err),
									 err.foc = length(err.foc),
									 err.sin = length(err.sin),
									 err.frq = paste(err.frq, collapse = ","))

	d = rbind(d, tmp)
}

tru.stack = d


if (file.exists(savefile)) next

save(est.stack, tru.stack, file = savefile)





