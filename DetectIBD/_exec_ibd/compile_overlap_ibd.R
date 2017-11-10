#
# compile overlapping IBD segments
#

library(data.table)


# args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" # args[1] # 
tru.file = "./result.truth.local.RData" # args[2] # 


load(sprintf("./data.%s.RData", prefix))

load(tru.file)



key = sprintf("%d %d %d %d", truth$h0, truth$h1, truth$lhs.index, truth$rhs.index)
del = which(duplicated(key))
tru = truth[-del, ]

key = sprintf("%d %d", tru$g0, tru$g1)
tab = table(key)
del = which(tab == 1)
del = names(tab)[del]
del = which(key %in% del)
tru = tru[-del, ]

key = sprintf("%d %d", tru$g0, tru$g1)
tru = tru[order(key), ]
key = sprintf("%d %d", tru$g0, tru$g1)
dup = duplicated(key)


lhs = tru$lhs.index[1]
rhs = tru$rhs.index[1]
rng = lhs:rhs

mat = rep(0, length(POS))
mat[rng] = 1

overlap = list()

k = 0
for (i in 2:nrow(tru)) {
	lhs = tru$lhs.index[i]
	rhs = tru$rhs.index[i]
	rng = lhs:rhs
	
	if (dup[i]) {
		
		tmp = rep(0, length(POS))
		tmp[rng] = 1
		mat = cbind(mat, tmp)
		
	} else {
		
		x = rowSums(mat)
		if (any(x > 1)) {
			overlap[[ key[i-1] ]] = mat
			
			if (length(overlap) == 1000) {
				cat("Saving ... \n")
				k = k + 1
				save(overlap, file = sprintf("result.overlap_ibd.%03d.RData", k))
				overlap = list()
				gc()
			}
		}
		
		mat = rep(0, length(POS))
		mat[rng] = 1
		
	}
	
	if (i %% 5000 == 0) cat(sprintf(" %d of %d, %d found\n", i, nrow(tru), length(overlap)))
}

k = k + 1
save(overlap, file = sprintf("result.overlap_ibd.%03d.RData", k))





###
stop("DONE")
###



library(data.table)


prefix = "OutOfAfricaHapMap20"

load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.G.RData", prefix)) # genotypes



cut.freq = function(f, brk = (0:500)/500, lab = sprintf("%.3f", (1:500)/500)) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

cut.geno = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}

make.matrix = function(len = nrow(G)) {
	matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
}

mask = c("00"=1, "01"=2, "02"=3, "11"=4, "12"=5, "22"=6)


res.files = dir(pattern = "^result\\.overlap_ibd\\..+\\.RData$")

count = make.matrix()

for (res.file in res.files) {
	cat(res.file)
	load(res.file)
	cat("\n")
	
	i = 0
	for (pair in names(overlap)) {
		i = i + 1
		if (i %% 100 == 0) cat(sprintf(" %d of %d\n", i, length(overlap)))
		
		g0 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", pair))
		g1 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", pair))
		
		g = cut.geno(G[, g0], G[, g1])
		
		x = rowSums(overlap[[pair]])
		x = which(x == 2)
		
		coord = matrix(c(x, mask[g[x]]), ncol = 2, byrow = F)
		count[ coord ] = count[ coord ] + 1
	}
}


save(count, file = sprintf("result.overlap_ibd_counts.%s.RData", prefix))






###
stop("DONE")
###



library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20"

load(sprintf("data.%s.RData", prefix))

load(sprintf("result.overlap_ibd_counts.%s.RData", prefix))


cut.freq = function(f, brk = (0:1000)/1000, lab = sprintf("%.3f", (1:1000)/1000)) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}


tmp = as.data.frame(count)
tmp = split(tmp, cut.freq(AAF))
tmp = lapply(tmp, function(x) {
	x = colSums(x)
	x = x / sum(x)
	as.data.table(as.list(x))
})
tag = names(tmp)
tmp = rbindlist(tmp)
rownames(tmp) = tag
freq = tmp



d = NULL

f = as.numeric(rownames(freq))

for (p in names(freq)) {
	x = freq[[p]]
	
	d = rbind(d, data.table(freq = f, 
													pair = p,
													prop = x))
}

del = which(is.na(d$prop))
if (length(del) > 0) {
	d = d[-del, ]
}




p = ((1:500)/500) - (0.5/500)
q = 1 - p

e = rbind(data.table(freq = q, pair = "00", prop = p^4),
					data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
					data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
					data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
					data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
					data.table(freq = q, pair = "22", prop = q^4))


expect.non = rbind(data.table(freq = q, pair = "00", prop = p^4),
									 data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
									 data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
									 data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
									 data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
									 data.table(freq = q, pair = "22", prop = q^4))
expect.non$type = "Expected under NON"

expect.ibd = rbind(data.table(freq = q, pair = "00", prop = p^3),
									 data.table(freq = q, pair = "01", prop = 2*(p^2)*q),
									 data.table(freq = q, pair = "02", prop = 0),
									 data.table(freq = q, pair = "11", prop = ((p^2)*q) + (p*(q^2))),
									 data.table(freq = q, pair = "12", prop = 2*p*(q^2)),
									 data.table(freq = q, pair = "22", prop = q^3))
expect.ibd$type = "Expected under IBD"



gg = ggplot(d) + 
	geom_point(aes(x = freq, y = prop, colour = pair), size=0.5) +
	geom_line(data = expect.ibd, aes(x = freq, y = prop, colour = pair), linetype = "22", alpha = 0.75) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("Allele frequency") +
	ylab("Observed genotype pair proportion") +
	ggtitle("Overlapping IBD")


ggsave(filename = "~/Desktop/genotypepairprops.overlap.pdf", plot = gg, width = 10, height = 10)




