

library(data.table)


infer.hap = function(tar, oth) {
	(tar + oth) / 4
}


load("../_exec_phasing/phase.rawdata.1043.RData")


## infer shared hap

infer.shared.hap = function(tar, oth, rng) {
	out = rep(NA, length(oth))
	hap = (tar[rng] + oth) / 4
	is0 = which(hap < 0.5)
	is1 = which(hap > 0.5)
	out[is0] = 0
	out[is1] = 1
	out
}

shap = list()

for (i in 1:nrow(coord)) {
	if (i %% 50 == 0) cat(".")
	if (i %% 5000 == 0) cat("\n")
	
	oth = chunk[[ coord$key[i] ]]
	rng = (coord$beg[i]):(coord$end[i])
	
	shap[[ coord$key[i] ]] = infer.shared.hap(target, oth, rng)
}
cat("\n")





focal = split(coord, coord$foc)


## fit to each other

fit.focal.hap = function(ha, hb, ra, rb, idx) {
	r = intersect(ra, rb)
	ia = match(r, ra)
	ib = match(r, rb)
	
	x = match(idx, r)
	
	d = abs(ha[ia] - hb[ib])
	
	dl = d[ 1:x ]
	dr = rev(d[ x:length(d) ])
	
	zl = min(r)
	zr = max(r)
	
	z = which(dl == 1); if (length(z) != 0) zl = zl + max(z)
	z = which(dr == 1); if (length(z) != 0) zr = zr - max(z)
	
	c(min(r), zl, max(r), zr)
}


k = 0
for (tag in names(focal)) {
	k = k + 1
	if (k %% 10 == 0) cat(".")
	if (k %% 1000 == 0) cat("\n")
	
	foc = focal[[tag]]
	
	foc$fit.beg = foc$beg
	foc$fit.end = foc$end
	
	idx = foc$foc[1]
	
	if (nrow(foc) != 1) {
	
	for (a in 1:nrow(foc)) {
		
		ha = shap[[ foc$key[a] ]]
		ra = (foc$beg[a]):(foc$end[a])
		
		fit = matrix(NA, nrow = nrow(foc), ncol = 4)
		
		for (b in 1:nrow(foc)) {
			if (a == b) next
			
			hb = shap[[ foc$key[b] ]]
			rb = (foc$beg[b]):(foc$end[b])
			
			fit[b, ] = fit.focal.hap(ha, hb, ra, rb, idx)
		}
		
		x = which(fit[, 1] == foc$beg[a])
		if (length(x) > 0) foc$fit.beg[a] = max(fit[x, 2], na.rm = T)
		
		x = which(fit[, 3] == foc$end[a])
		if (length(x) > 0) foc$fit.end[a] = min(fit[x, 4], na.rm = T)
	}
	}
	
	focal[[tag]] = foc
}
cat("\n")

fitted = lapply(focal, as.data.table)
fitted = rbindlist(fitted)




#
# stack
#

sorted = fitted[order(fitted$fit.end), ]

sorted$infhap = NA
sorted$infhap[1] = sorted$h0[1]

for (i in 2:nrow(sorted)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")
	
	if (sorted$h0[i] == sorted$h1[i]) next
	
	prevhap = sorted$infhap[i-1]
	
	beg = seg$fit.beg[i]
	end = seg$fit.end[i]
	
	ovr = which(sorted$fit.beg <= end & end <= sorted$fit.end)
	if (length(ovr) == 0) next
	ovr = fitted[ovr, ]
	
	
}




#
# distinguish haplotypes
#

het = which(target == 1)

distinguish.hap = function(ha, hb, ra, rb, het) {
	r = intersect(ra, rb)
	l = length(r)
	r = intersect(r, het)
	ia = match(r, ra)
	ib = match(r, rb)
	
	d = abs(ha[ia] - hb[ib])
	d = table(d, useNA = "always")
	
	is0 = d["0"]; if(is.na(is0)) is0 = 0
	is1 = d["1"]; if(is.na(is1)) is1 = 0
	isX = sum(d) - (is0 + is1)
	
	c(l, length(r), isX, is0 + is1, is0, is1)
}


seg = sample(1:nrow(fitted), min(nrow(fitted), 100))
seg = fitted[seg, ]


# before fitting

d = NULL

for (i in 1:nrow(seg)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")
	
	beg = seg$beg[i]
	end = seg$end[i]
	
	if (seg$h0[i] == seg$h1[i]) next
	
	ovr = which((fitted$beg <= beg & beg <= fitted$end) | (fitted$beg <= end & end <= fitted$end))
	if (length(ovr) == 0) next
	ovr = fitted[ovr, ]
	
	ref = shap[[ seg$key[i] ]]
	rr = beg:end
	rh = seg$h0[i]
	
	os = which(ovr$h0)
	od = which(!ovr$h0)
	if (length(os) > 50) os = sample(os, 50)
	if (length(od) > 50) od = sample(od, 50)
	ovr = ovr[c(os, od), ]
	
	dis = data.table(fk.seg = seg$fk[i], fk.ovr = ovr$fk, hap = (rh == ovr$h0))
	tmp = matrix(0, nrow = nrow(ovr), ncol = 6, dimnames = list(NULL, c("len", "het", "hht", "inf", "same", "diff")))
	
	for (j in 1:nrow(ovr)) {
		if (ovr$h0[j] == ovr$h1[j]) next
		
		alt = shap[[ ovr$key[j] ]]
		aa = (ovr$beg[j]):(ovr$end[j])
		ah = ovr$h0[j]
		
		tmp[j, ] = distinguish.hap(ref, alt, rr, aa, het)
	}
	
	tmp = as.data.table(tmp)
	dis = cbind(dis, tmp)
	
	del = which(is.na(dis$len))
	if (length(del) > 0) {
		dis = dis[-del, ]
	}
	
	d = rbind(d, dis)
}
cat("\n")

est.dis = d



# after fitting

d = NULL

for (i in 1:nrow(seg)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")
	
	beg = seg$fit.beg[i]
	end = seg$fit.end[i]
	
	if (seg$h0[i] == seg$h1[i]) next
	
	ovr = which((fitted$fit.beg <= beg & beg <= fitted$fit.end) | (fitted$fit.beg <= end & end <= fitted$fit.end))
	if (length(ovr) == 0) next
	ovr = fitted[ovr, ]
	
	ref = shap[[ seg$key[i] ]]
	zz = (seg$beg[i]):(seg$end[i])
	rr = intersect(beg:end, zz)
	if (length(rr) == 0) stop("???")
	ref = ref[match(rr, zz)]
	rh = seg$h0[i]
	
	os = which(ovr$h0)
	od = which(!ovr$h0)
	if (length(os) > 50) os = sample(os, 50)
	if (length(od) > 50) od = sample(od, 50)
	ovr = ovr[c(os, od), ]
	
	dis = data.table(fk.seg = seg$fk[i], fk.ovr = ovr$fk, hap = (rh == ovr$h0))
	tmp = matrix(0, nrow = nrow(ovr), ncol = 6, dimnames = list(NULL, c("len", "het", "hht", "inf", "same", "diff")))
	
	for (j in 1:nrow(ovr)) {
		if (ovr$h0[j] == ovr$h1[j]) next
		
		alt = shap[[ ovr$key[j] ]]
		zz = (ovr$beg[j]):(ovr$end[j])
		aa = intersect((ovr$fit.beg[j]):(ovr$fit.end[j]), zz)
		if (length(aa) == 0) stop("??????")
		alt = alt[match(aa, zz)]
		ah = ovr$h0[j]
		
		tmp[j, ] = distinguish.hap(ref, alt, rr, aa, het)
	}
	
	tmp = as.data.table(tmp)
	dis = cbind(dis, tmp)
	
	del = which(is.na(dis$len))
	if (length(del) > 0) {
		dis = dis[-del, ]
	}
	
	d = rbind(d, dis)
}
cat("\n")

fit.dis = d



# true segments

d = NULL

for (i in 1:nrow(seg)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")
	
	beg = seg$true.beg[i]
	end = seg$true.end[i]
	
	if (seg$h0[i] == seg$h1[i]) next
	
	ovr = which((fitted$true.beg <= beg & beg <= fitted$true.end) | (fitted$true.beg <= end & end <= fitted$true.end))
	if (length(ovr) == 0) next
	ovr = fitted[ovr, ]
	
	ref = shap[[ seg$key[i] ]]
	zz = (seg$beg[i]):(seg$end[i])
	rr = intersect(beg:end, zz)
	if (length(rr) == 0) stop("???")
	ref = ref[match(rr, zz)]
	rh = seg$h0[i]
	
	os = which(ovr$h0)
	od = which(!ovr$h0)
	if (length(os) > 50) os = sample(os, 50)
	if (length(od) > 50) od = sample(od, 50)
	ovr = ovr[c(os, od), ]
	
	dis = data.table(fk.seg = seg$fk[i], fk.ovr = ovr$fk, hap = (rh == ovr$h0))
	tmp = matrix(0, nrow = nrow(ovr), ncol = 6, dimnames = list(NULL, c("len", "het", "hht", "inf", "same", "diff")))
	
	for (j in 1:nrow(ovr)) {
		if (ovr$h0[j] == ovr$h1[j]) next
		
		alt = shap[[ ovr$key[j] ]]
		zz = (ovr$beg[j]):(ovr$end[j])
		aa = intersect((ovr$true.beg[j]):(ovr$true.end[j]), zz)
		if (length(aa) == 0) stop("??????")
		alt = alt[match(aa, zz)]
		ah = ovr$h0[j]
		
		tmp[j, ] = distinguish.hap(ref, alt, rr, aa, het)
	}
	
	tmp = as.data.table(tmp)
	dis = cbind(dis, tmp)
	
	del = which(is.na(dis$len))
	if (length(del) > 0) {
		dis = dis[-del, ]
	}
	
	d = rbind(d, dis)
}
cat("\n")

tru.dis = d




###
### clustering
###

het = which(target == 1)

M = matrix(NA, nrow(coord), length(het))

for (i in 1:nrow(coord)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")
	
	beg = coord$beg[i]
	end = coord$end[i]
	
	rng = beg:end
	rrr = match(het, rng)
	
	shp = shap[[ coord$key[i] ]]
	
	M[i, ] = shp[rrr]
}

rownames(M) = coord$key

M = M[sample(coord$key, 500), ]
del = which(apply(M, 2, function(x) all(is.na(x))))
if (length(del) > 0) {
	M = M[, -del]
}

D = dist(M, method = "euclidean")
C = hclust(D)
K = cutree(C, k = 2)





stop()
##############

est.dis$mode = "Before"
fit.dis$mode = "After"

d = rbind(est.dis, fit.dis)
d$mode = factor(d$mode, levels = c("Before", "After"), ordered = T)

d$prop  = log10(d$inf)
d$delta = (d$same - d$diff)
d$class = "Indistinguishable"

x = which(d$delta < 0 & !d$hap); d$class[x] = "Correct"
x = which(d$delta < 0 & d$hap);  d$class[x] = "Incorrect"

x = which(d$delta > 0 & d$hap);  d$class[x] = "Correct"
x = which(d$delta > 0 & !d$hap); d$class[x] = "Incorrect"

x = which(d$delta < 0); d$delta[x] = log10(abs(d$delta[x])) * -1
x = which(d$delta > 0); d$delta[x] = log10(d$delta[x])

d$class = factor(d$class, levels = c("Correct", "Incorrect", "Indistinguishable"), ordered = T)
d = d[order(d$class), ]


ggplot(d) + 
	facet_grid(.~mode) + 
	geom_point(aes(x=prop, y=delta, colour = hap, shape = class)) + 
	scale_shape_manual(values = c("Correct" = 1, "Incorrect" = 16, "Indistinguishable" = 4))


ggplot(d) + 
	facet_grid(.~mode) + 
	geom_bin2d(aes(x=prop, y=delta), binwidth = c(1/30, 1/30)) +
	theme(aspect.ratio = 2)


ggplot(d) +
	facet_grid(.~mode) + 
	geom_bar(aes(x = fk.ovr, fill = class), position = "fill")







