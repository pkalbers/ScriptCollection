

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
