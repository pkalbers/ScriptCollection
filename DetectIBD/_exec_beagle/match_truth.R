#
# include Beagle IBD results
#

library(data.table)


prefix = "OutOfAfricaHapMap20"


ibd.H.file = sprintf("beagle_ibd.%s.H.ibd", prefix)
ibd.P.file = sprintf("beagle_ibd.%s.P.ibd", prefix)

tru.file = "./result.truth.local.RData"


load(sprintf("data.%s.RData", prefix))
load(tru.file)


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



ibd.H = read.table(ibd.H.file, header = F, stringsAsFactors = F)
ibd.P = read.table(ibd.P.file, header = F, stringsAsFactors = F)

names(ibd.H) = c("indv1", "x1", "indv2", "x2", "chr", "beg", "end", "lod")
names(ibd.P) = c("indv1", "x1", "indv2", "x2", "chr", "beg", "end", "lod")



ibd.H$g0 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd.H$indv1))
ibd.H$g1 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd.H$indv2))

ibd.P$g0 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd.P$indv1))
ibd.P$g1 = as.numeric(sub("^ID([0-9]+)$", "\\1", ibd.P$indv2))


ibd.H$h0 = (ibd.H$g0 * 2) + (ibd.H$x1 - 2)
ibd.H$h1 = (ibd.H$g1 * 2) + (ibd.H$x2 - 2)

ibd.P$h0 = (ibd.P$g0 * 2) + (ibd.P$x1 - 2)
ibd.P$h1 = (ibd.P$g1 * 2) + (ibd.P$x2 - 2)


beg = match(ibd.H$beg, pos)
end = match(ibd.H$end, pos)

ibd.H$beg.pos = POS[beg]
ibd.H$end.pos = POS[end]

ibd.H$beg.idx = beg
ibd.H$end.idx = end


beg = match(ibd.P$beg, pos)
end = match(ibd.P$end, pos)

ibd.P$beg.pos = POS[beg]
ibd.P$end.pos = POS[end]

ibd.P$beg.idx = beg
ibd.P$end.idx = end



ibd = ibd.H  ### <-- MANUALLY CHANGE HERE
tru = truth

tru.key = sprintf("%d %d", tru$h0, tru$h1)
ibd.key = sprintf("%d %d", ibd$h0, ibd$h1)

key = intersect(ibd.key, tru.key)

tru = tru[which(tru.key %in% key), ]
ibd = ibd[which(ibd.key %in% key), ]

tru.key = sprintf("%d %d", tru$h0, tru$h1)
ibd.key = sprintf("%d %d", ibd$h0, ibd$h1)

match = NULL

while(T) {
	key = intersect(ibd.key, tru.key)
	
	if (length(key) == 0) {
		break
	}
	
	tru.idx = match(key, tru.key)
	ibd.idx = match(key, ibd.key)
	
	sub.tru = tru[tru.idx, ]
	sub.ibd = ibd[ibd.idx, ]
	
	a = which(sub.tru$position >= sub.ibd$beg.pos)
	b = which(sub.tru$position <  sub.ibd$end.pos)
	x = intersect(a, b)
	
	if (length(x) == 0) {
		if (length(ibd.idx) == nrow(ibd)) {
			break
		}
		cat(sprintf("Reset to %d\n", nrow(ibd)))
		ibd = ibd[-(ibd.idx), ]
		ibd.key = ibd.key[-(ibd.idx)]
		next
	}
	
	cat(sprintf(" %d of %d\n", length(x), nrow(tru)))
	
	tmp = as.data.table(sub.tru[x, ])
	tmp = cbind(tmp, data.table(beagle.h0 = sub.ibd$h0[x],
															beagle.h1 = sub.ibd$h1[x],
															beagle.lhs.position = sub.ibd$beg.pos[x],
															beagle.rhs.position = sub.ibd$end.pos[x],
															beagle.lhs.index = sub.ibd$beg.idx[x],
															beagle.rhs.index = sub.ibd$end.idx[x]))
	match = rbind(match, tmp)
	
	tru = tru[-(tru.idx[x]), ]
	tru.key = tru.key[-(tru.idx[x])]
}


print(all(match$h0 == match$beagle.h0))
print(all(match$h1 == match$beagle.h1))

match.H = match  ### <-- MANUALLY CHANGE HERE




ibd = ibd.P  ### <-- MANUALLY CHANGE HERE
tru = truth

tru.key = sprintf("%d %d", tru$h0, tru$h1)
ibd.key = sprintf("%d %d", ibd$h0, ibd$h1)

key = intersect(ibd.key, tru.key)

tru = tru[which(tru.key %in% key), ]
ibd = ibd[which(ibd.key %in% key), ]

tru.key = sprintf("%d %d", tru$h0, tru$h1)
ibd.key = sprintf("%d %d", ibd$h0, ibd$h1)

match = NULL

while(T) {
	key = intersect(ibd.key, tru.key)
	
	if (length(key) == 0) {
		break
	}
	
	tru.idx = match(key, tru.key)
	ibd.idx = match(key, ibd.key)
	
	sub.tru = tru[tru.idx, ]
	sub.ibd = ibd[ibd.idx, ]
	
	a = which(sub.tru$position >= sub.ibd$beg.pos)
	b = which(sub.tru$position <  sub.ibd$end.pos)
	x = intersect(a, b)
	
	if (length(x) == 0) {
		if (length(ibd.idx) == nrow(ibd)) {
			break
		}
		cat(sprintf("Reset to %d\n", nrow(ibd)))
		ibd = ibd[-(ibd.idx), ]
		ibd.key = ibd.key[-(ibd.idx)]
		next
	}
	
	cat(sprintf(" %d of %d\n", length(x), nrow(tru)))
	
	tmp = as.data.table(sub.tru[x, ])
	tmp = cbind(tmp, data.table(beagle.h0 = sub.ibd$h0[x],
															beagle.h1 = sub.ibd$h1[x],
															beagle.lhs.position = sub.ibd$beg.pos[x],
															beagle.rhs.position = sub.ibd$end.pos[x],
															beagle.lhs.index = sub.ibd$beg.idx[x],
															beagle.rhs.index = sub.ibd$end.idx[x]))
	match = rbind(match, tmp)
	
	tru = tru[-(tru.idx[x]), ]
	tru.key = tru.key[-(tru.idx[x])]
}


print(all(match$h0 == match$beagle.h0))
print(all(match$h1 == match$beagle.h1))

match.P = match  ### <-- MANUALLY CHANGE HERE




save(match.H, match.P, file = sprintf("result.beagle.%s.RData", prefix))









