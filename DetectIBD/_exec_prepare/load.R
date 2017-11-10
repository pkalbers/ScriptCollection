#
# load simulated data into R
#

library(data.table)


prefix = "OutOfAfricaHapMap20"

samplesize = 5000


# get simulated haplotypes
cat("Reading haplotypes\n")

con = file(sprintf("%s.hdf5.hap", prefix), "r", blocking = FALSE)
line = readLines(con, n = 1)
nmarker = nchar(line)
close(con)

H = matrix(integer(), nrow = nmarker, ncol = samplesize)

con = file(sprintf("%s.hdf5.hap", prefix), "r", blocking = FALSE)
nsample = 0
while(length(line <- readLines(con, n = 1)) == 1) {
	nsample = nsample + 1
	line = strsplit(line, "", T)[[1]]
	line = as.integer(line)
	H[, nsample] = line
	if (nsample %% 100 == 0) cat(sprintf(" %d of %d\n", nsample, samplesize))
}
cat("\n")
close(con)

cat("# haplotypes: ", ncol(H), "\n")
cat("# markers:    ", nrow(H), "\n")


save(H, file = sprintf("data.%s.H.RData", prefix))


# make genotypes
cat("Making genotypes\n")
G = matrix(integer(), nrow = nrow(H), ncol = ncol(H) / 2)
j = 1
for (i in 1:ncol(G)) {
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
	if (i %% 50 == 0) cat(sprintf(" %d of %d\n", i, ncol(G)))
}
cat("\n")


save(G, file = sprintf("data.%s.G.RData", prefix))


# get frequencies
cat("Getting frequencies\n")
MAC = AAC = rowSums(H)
MAF = AAF = AAC / ncol(H)
i = which(MAC > (ncol(H) / 2))
if (length(i) > 0) {
	MAC[i] = ncol(H) - MAC[i]
	MAF[i] = 1 - MAF[i]
}

GC0 = apply(G, 1, function(x) length(which(x == 0)))
GC1 = apply(G, 1, function(x) length(which(x == 1)))
GC2 = apply(G, 1, function(x) length(which(x == 2)))

GF0 = GC0 / ncol(G)
GF1 = GC1 / ncol(G)
GF2 = GC2 / ncol(G)


# get positions
POS = as.numeric(readLines(sprintf("%s.hdf5.pos", prefix)))
CPOS = readLines(sprintf("%s.hdf5.pos", prefix))

# get breakpoints
BRK = as.numeric(readLines(sprintf("%s.hdf5.brk", prefix)))
CBRK = readLines(sprintf("%s.hdf5.brk", prefix))


# make fk index
cat("Compiling rare variant index\n")
fklist = 2:25 # c(2:25, seq(30, 45, by = 5), seq(50, 90, by = 10), seq(100, 500, by=50))
FKI = list()
cum = 0
i = 0
for (fk in fklist) {
	i = i + 1
	idx = which(AAC == fk)
	num = length(idx)
	cmb = choose(fk, 2) * num
	cum = cmb + cum
	cat(sprintf("# f_%d : %d (%d pairs | %d cumulative pairs)\n", fk, num, cmb, cum))
	
	tmp = lapply(idx, function(x) {
		g = which(G[x, ] == 1)
		if (length(g) < 2) return(NULL) # others may be autozygous, hence higher count
		as.data.table(cbind(x, POS[x], AAC[x], t(combn(g, 2))))
	})
	tmp = rbindlist(tmp)
	names(tmp) = c("index", "pos", "fk", "g0", "g1")
	
	a = tmp$g0 * 2 - 1
	b = tmp$g0 * 2
	coord.a = as.matrix(data.frame(x = tmp$index, y = a))
	coord.b = as.matrix(data.frame(x = tmp$index, y = b))
	coord.a = H[ coord.a ]
	coord.b = H[ coord.b ]
	if (any(coord.a == coord.b)) stop("???")
	tmp$h0 = (a * coord.a) + (b * coord.b)
	
	a = tmp$g1 * 2 - 1
	b = tmp$g1 * 2
	coord.a = as.matrix(data.frame(x = tmp$index, y = a))
	coord.b = as.matrix(data.frame(x = tmp$index, y = b))
	coord.a = H[ coord.a ]
	coord.b = H[ coord.b ]
	if (any(coord.a == coord.b)) stop("???")
	tmp$h1 = (a * coord.a) + (b * coord.b)
	
	FKI = rbind(FKI, tmp)
}


# fklist = c(2:(samplesize - 2))
# FKD = data.table(fk = fklist, count = NA, pairs = NA, cum.pairs = NA)
# cum = 0
# i = 0
# for (fk in fklist) {
# 	i = i + 1
# 	idx = which(AAC == fk)
# 	num = length(idx)
# 	cmb = choose(fk, 2) * num
# 	cum = cmb + cum
# 
# 	FKD$count[i] = num
# 	FKD$pairs[i] = cmb
# 	FKD$cum.pairs[i] = cum
# }


# save data
cat("Saving data\n")
save(POS, CPOS, BRK, CBRK, 
		 AAC, AAF, MAC, MAF, 
		 GC0, GC1, GC2, GF0, GF1, GF2,
		 FKI, # FKD,
		 file = sprintf("data.%s.RData", prefix))

cat("OK\n")


#
# write sharing pair table for IBD detection
#

cat("Distributing sharing pair tables ... \n")

dir.create("./_truth", showWarnings = F)

tmp = split(FKI, cut(1:nrow(FKI), 500))

for (i in 1:length(tmp)) {
	ibd.name = sprintf("./_truth/pairs.%04d.%s.txt", i, prefix)
	ibd.file = file(ibd.name, open = "w")
	pair = tmp[[i]]
	for (x in 1:nrow(pair)) {
		cat(sprintf("%d %s %d %d %d %d %d\n", pair$index[x], CPOS[pair$index[x]], pair$fk[x], pair$g0[x], pair$g1[x], pair$h0[x], pair$h1[x]), file = ibd.file)
	}
	close(ibd.file)
}

cat("\nDONE!\n")












