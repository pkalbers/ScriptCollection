

infer.hap = function(tar, oth) {
	(tar + oth) / 4
}


files = dir(pattern = "^phase\\.rawdata\\.(.+)\\.RData$", path = "../rawdata")
files = rep(files, 10)

for (file in sample(files)) {
	x = sub("^phase\\.rawdata\\.(.+)\\.RData$", "\\1", file)

	TE = sample(c("est", "tru"), 1)
	
	flagfile = sprintf("_flag.overlap.%s.%s.txt", x, TE)
	savefile = sprintf("phase.overlap.%s.%s.RData", x, TE)

	if (file.exists(savefile)) next

	Sys.sleep(runif(1, min = 0.01, max = 1))

	if (file.exists(flagfile)) next

	cat(".", file = flagfile, append = F)

	load(sprintf("../rawdata/%s", file))

	if (idv != as.numeric(x)) stop(sprintf("%d %s", idv, x))

	cat("\n\n", idv, "\n\n")


	coord = cbind(this = idv, coord)


#
# overlap, estimated length
#
if (TE == "est") {
	
est.covered = T

# test complete cover
cov = rep(F, max(max(coord$beg), max(coord$true.beg)))
for (i in 1:nrow(coord)) {
	if (i %% 100 == 0) cat(".")
	if (i %% 5000 == 0) cat("\n")

	a = coord$beg[i]
	b = coord$end[i]
	rng = a:b

	cov[rng] = T
}
cat("\n")

if (any(!cov)) {
	est.covered = F
}


over = coord
over = over[order(over$end), ]
rownames(over) = as.character(1:nrow(over))
over = split(over, 1:nrow(over))

overlap = lapply(over, function(over, coord, chunk) {
	cat(".")
	
	foc = over$foc
	beg = over$beg
	end = over$end
	
	ovr = which(coord$foc != foc & coord$beg < end & end < coord$end)
	
	if (length(ovr) == 0) {
		return(NULL)
	}
	
	ref = chunk[[ over$key ]]
	ref.rng = beg:end
	
	sub = coord[ovr, ]
	sub$ref = as.numeric(rownames(over))
	sub$hap = NA
	sub$len = 0
	sub$het = 0
	sub$inf.ref = 0
	sub$inf.alt = 0
	sub$inf.int = 0
	sub$dis = 0
	sub$dis.frq = ""
	
	for (j in 1:nrow(sub)) {
		if (over$h0 == sub$h0[j]) {
			sub$hap[j] = T
		} else {
			sub$hap[j] = F
		}
		
		alt = chunk[[ sub$key[j] ]]
		alt.rng = (sub$beg[j]):(sub$end[j])
		
		abs.rng = intersect(ref.rng, alt.rng)
		
		if (length(abs.rng) == 0) next
		
		tar = target[abs.rng]
		het = which(tar == 1)
		
		rr = match(abs.rng, ref.rng)
		aa = match(abs.rng, alt.rng)
		
		hr = infer.hap(tar, ref[rr])
		ha = infer.hap(tar, alt[aa])
		
		inf.ref = which(hr[het] != 0.5)
		inf.alt = which(ha[het] != 0.5)
		inf.int = intersect(inf.ref, inf.alt)
		
		dis = which(abs(hr - ha) >= 0.5)
		
		dis.frq = c()
		if (length(dis) > 0) {
			dis.frq = frq[abs.rng[dis]]
		}
		
		sub$len[j] = length(abs.rng)
		sub$het[j] = length(het)
		sub$inf.ref[j] = length(inf.ref)
		sub$inf.alt[j] = length(inf.alt)
		sub$inf.int[j] = length(inf.int)
		sub$dis[j] = length(dis)
		sub$dis.frq[j] = paste(dis.frq, collapse = ",")
	}
	
	sub$key = NULL
	sub
}, coord, chunk)
cat("\n")

est.overlap = overlap

save(est.covered, est.overlap, file = savefile)
}

#
# overlap, true length
#
if (TE == "tru") {
tru.covered = T

# test complete cover
cov = rep(F, max(max(coord$beg), max(coord$true.beg)))
for (i in 1:nrow(coord)) {
	if (i %% 100 == 0) cat(".")
	if (i %% 5000 == 0) cat("\n")

	a = coord$true.beg[i]
	b = coord$true.end[i]
	rng = a:b

	cov[rng] = T
}
cat("\n")

if (any(!cov)) {
	tru.covered = F
}


over = coord
over = over[order(over$true.end), ]
rownames(over) = as.character(1:nrow(over))
over = split(over, 1:nrow(over))

overlap = lapply(over, function(over, coord, chunk) {
	cat(".")
	
	foc = over$foc
	beg = over$true.beg
	end = over$true.end
	
	if (beg < over$beg || end > over$end) return(NULL)
	
	ovr = which(coord$foc != foc & coord$true.beg < end & end < coord$true.end)
	
	if (length(ovr) == 0) {
		return(NULL)
	}
	
	ref = chunk[[ over$key ]]
	ref.rng = beg:end
	
	tmp = (over$beg):(over$end)
	del = which(is.na(match(tmp, ref.rng)))
	if (length(del) > 0) ref = ref[-del]
	
	sub = coord[ovr, ]
	sub$ref = as.numeric(rownames(over))
	sub$hap = NA
	sub$len = 0
	sub$het = 0
	sub$inf.ref = 0
	sub$inf.alt = 0
	sub$inf.int = 0
	sub$dis = 0
	sub$dis.frq = ""
	
	for (j in 1:nrow(sub)) {
		if (over$h0 == sub$h0[j]) {
			sub$hap[j] = T
		} else {
			sub$hap[j] = F
		}
		
		alt = chunk[[ sub$key[j] ]]
		alt.rng = (sub$true.beg[j]):(sub$true.end[j])
		
		tmp = (sub$beg[j]):(sub$end[j])
		del = which(is.na(match(tmp, alt.rng)))
		if (length(del) > 0) alt = alt[-del]
		
		abs.rng = intersect(ref.rng, alt.rng)
		
		if (length(abs.rng) == 0) next
		
		tar = target[abs.rng]
		het = which(tar == 1)
		
		rr = match(abs.rng, ref.rng)
		aa = match(abs.rng, alt.rng)
		
		hr = infer.hap(tar, ref[rr])
		ha = infer.hap(tar, alt[aa])
		
		inf.ref = which(hr[het] != 0.5)
		inf.alt = which(ha[het] != 0.5)
		inf.int = intersect(inf.ref, inf.alt)
		
		dis = which(abs(hr - ha) >= 0.5)
		
		dis.frq = c()
		if (length(dis) > 0) {
			dis.frq = frq[abs.rng[dis]]
		}
		
		sub$len[j] = length(abs.rng)
		sub$het[j] = length(het)
		sub$inf.ref[j] = length(inf.ref)
		sub$inf.alt[j] = length(inf.alt)
		sub$inf.int[j] = length(inf.int)
		sub$dis[j] = length(dis)
		sub$dis.frq[j] = paste(dis.frq, collapse = ",")
	}
	
	sub$key = NULL
	sub
}, coord, chunk)
cat("\n")

tru.overlap = overlap

save(tru.covered, tru.overlap, file = savefile)
}
}


stop()
###

d = lapply(overlap, function(x) {
	del = which(x$inf == 0)
	if (length(del) > 0) {
		x = x[-del,]
	}
	if (nrow(x) == 0) return(NULL)
	sapply(split(x, x$hap), function(y) {
		y$dis / y$inf
	})
})
dt = unlist(lapply(d, function(x) if (!"TRUE" %in% names(x)) return(NULL) else x[["TRUE"]]), use.names = F)
df = unlist(lapply(d, function(x) if (!"FALSE" %in% names(x)) return(NULL) else x[["FALSE"]]), use.names = F)






