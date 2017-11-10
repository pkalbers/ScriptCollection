


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]


load("result.match.OutOfAfricaHapMap20.RData")


load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.H.RData", prefix)) # original haplotypes
load(sprintf("data.%s.G.RData", prefix)) # original genotypes


IDX = 1:length(POS)
names(IDX) = as.character(POS)




#
# sharing
#
if (FALSE) {
	share = rep(0, ncol(G))
	for (i in 1:ncol(G)) {
		if (i %% 10 == 0) cat(".")
		if (i %% 1000 == 0) cat("\n")
		share[i] = length(which(match$g0 == i)) + length(which(match$g1 == i))
	}

	se <- function(x) sqrt(var(x)/length(x))

	mean(share)
	se(share)
}




# singletons
frq = rowSums(H)
singletons = which(frq == 1)




### target
Sys.sleep(runif(n = 1, min = 1, max = 5))
qwerty = 0
while (qwerty < 500) {
	qwerty = qwerty + 1


idv = sample(1:ncol(G), 1)


savename = sprintf("./rawdata/phase.rawdata.%04d.RData", idv)
flagname = sprintf("./rawdata/_flag.phase.rawdata.%04d.txt", idv)


if (file.exists(savename)) next

Sys.sleep(runif(n = 1, min = 0.01, max = 1))

if (file.exists(flagname)) next

cat(".", file = flagname, append = F)



m = which(match$g0 == idv | match$g1 == idv)
m = match[m, ]

target = G[, idv]
target.h0 = H[, (idv*2 - 1)]
target.h1 = H[, (idv*2)]

#target.het = which(target == 1)


# singletons
single = intersect(which(target == 1), singletons)



#
# coverage
#
if (file.exists(savename)) next

cov = rep(0, length(POS))
cov.fk = lapply(sort(unique(m$fk)), function(x) rep(0, length(POS)))
names(cov.fk) = as.character(sort(unique(m$fk)))

for (i in 1:nrow(m)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")

	a = IDX[ as.character(m$g.lhs[i]) ]
	b = IDX[ as.character(m$g.rhs[i]) ]

	if (any(is.na(a)) || any(is.na(b))) stop("???")

	rng = a:b

	cov[rng] = cov[rng] + 1

	fk = as.character(m$fk[i])
	cov.fk[[ fk ]][rng] = cov.fk[[ fk ]][rng] + 1
}
cat("\n")

range(cov)


# cumulative
cov.cum = list()

for (fk in sort(unique(m$fk))) {
	tag = as.character(fk)
	cum = rep(0, length(POS))

	for (sub in 2:fk) {
		cum = cum + cov.fk[[ as.character(sub) ]]
	}

	cov.cum[[tag]] = cum
}



d = NULL # data.frame(fk = "all", idx = IDX, cov = cov)

for (tag in names(cov.cum)) {
	cat(tag, " ")
	tmp = cov.cum[[tag]]

	if (length(tmp) == 0) next

	del = rep(F, length(POS))
	for (i in 2:(length(POS) - 1)) {
		if (tmp[i-1] == tmp[i] && tmp[i] == tmp[i+1]) del[i] = T
	}
	del = which(del)

	d = rbind(d, data.frame(fk = as.numeric(tag), idx = IDX[-del], pos = POS[-del], cov = tmp[-del]))
}
cat("\n")

total.cover = cov
cover = d

# d$cov = log10(d$cov)
# if (any(d$cov < 0)) d$cov[which(d$cov < 0)] = 0
#
# ggplot(d) + geom_step(aes(idx, cov, colour = factor(fk)))



#
# pair count
#
if (file.exists(savename)) next

num = matrix(0, nrow = nrow(m), ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))

for (i in 1:nrow(m)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")

	oth = m$g0[i]
	if (idv == oth) oth = m$g1[i]

	other = G[, oth]

	a = IDX[ as.character(m$g.lhs[i]) ]
	b = IDX[ as.character(m$g.rhs[i]) ]

	rng = a:b

	tar = target[rng]
	oth = other[rng]

	pair = sprintf("%d%d", tar, oth)
	x = which(pair == "10"); if (length(x) > 0) pair[x] = "01"
	x = which(pair == "20"); if (length(x) > 0) pair[x] = "02"
	x = which(pair == "21"); if (length(x) > 0) pair[x] = "12"

	pair = factor(pair, levels = c("00", "01", "02", "11", "12", "22"))

	num[i, ] = as.vector(table(pair))
}
cat("\n")

count = num



#
# collect sequence chunks
#
if (file.exists(savename)) next

chunk.crd = data.frame(other = 0, fk = 0, foc = 0, key = "", beg = rep(0, nrow(m)), end = rep(0, nrow(m)), true.beg = rep(0, nrow(m)), true.end = rep(0, nrow(m)), stringsAsFactors = F)
chunk.seq = list()

for (i in 1:nrow(m)) {
	if (i %% 10 == 0) cat(".")
	if (i %% 1000 == 0) cat("\n")

	key = sprintf("%d-%d-%d", m$index[i], m$g0[i], m$g1[i])

	oth = m$g0[i]
	chunk.crd$other[i] = oth
	if (idv == oth) {
		oth = m$g1[i]
		chunk.crd$other[i] = oth
	}

	other = G[, oth]

	a = IDX[ as.character(m$g.lhs[i]) ]
	b = IDX[ as.character(m$g.rhs[i]) ]

	rng = a:b

	chunk = other[rng]

	chunk.crd$fk[i] = m$fk[i]
	chunk.crd$foc[i] = m$index[i]
	chunk.crd$key[i] = key
	chunk.crd$beg[i] = a
	chunk.crd$end[i] = b
	chunk.crd$true.beg[i] = m$true.lhs.idx[i]
	chunk.crd$true.end[i] = m$true.rhs.idx[i]

	chunk.seq[[key]] = chunk
}
cat("\n")

chunk.crd$h0 = (target.h0[ chunk.crd$foc ] == 1)
chunk.crd$h1 = (target.h1[ chunk.crd$foc ] == 1)

coord = chunk.crd
chunk = chunk.seq


if (file.exists(savename)) next

save(idv, frq, single, target, target.h0, target.h1, total.cover, cover, count, coord, chunk, file = savename)

}
















