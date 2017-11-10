

load("../frq.RData")


flip.hap = function(h) {
	x = which(!is.na(h))
	h[x] = abs(h[x] - 1)
	h
}


infer.hap = function(tar, oth) {
	(tar + oth) / 4
}


files = dir(pattern = "^phase\\.shap\\.(.+)\\.RData$", path = "../shap")

for (file in sample(files)) {
	x = sub("^phase\\.shap\\.(.+)\\.RData$", "\\1", file)

	flagfile = sprintf("_flag.stack.%s.txt", x)
	savefile = sprintf("phase.stack.%s.RData", x)

	if (file.exists(savefile)) next

	Sys.sleep(runif(1, min = 0.01, max = 1))

	if (file.exists(flagfile)) next

	cat(".", file = flagfile, append = F)

	load(sprintf("../shap/%s", file))



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


}



stop()
###

library(data.table)

files = dir(pattern = "^phase\\.stack\\.(.+)\\.RData$")

es = NULL
ts = NULL

for (file in files) {
	cat(".")
	load(file)
	
	est.stack$err.frq = as.character(est.stack$err.frq)
	tru.stack$err.frq = as.character(tru.stack$err.frq)
	
	est.stack = as.data.table(est.stack)
	tru.stack = as.data.table(tru.stack)
	
	es = rbind(es, est.stack)
	ts = rbind(ts, tru.stack)
}
cat("\n")

est.stack = es
tru.stack = ts

save(est.stack, tru.stack, file = "../phase.stack.all.RData")











