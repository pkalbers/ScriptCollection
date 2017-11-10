#
# load simulated data into R
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

fk.max = as.numeric(args[1]) # rare variant max count
prefix = args[2]             # file prefix
shapairs = as.numeric(args[3]) # number of distributed sharing pairs files



# get simulated haplotypes
cat("Reading haplotypes\n")

con = file(sprintf("%s.VarByCol.hap", prefix), "r", blocking = FALSE)
nmarker = 0
nsample = 0
while(length(line <- readLines(con, n = 1)) == 1) {
	if (nmarker == 0) {
		nmarker = nchar(line)
	} else {
		if (nmarker != nchar(line)) stop("Different number of markers")
	}
	nsample = nsample + 1
	if (nsample %% 100 == 0) cat(".")
}
cat("\n")
close(con)

H = matrix(integer(), nrow = nmarker, ncol = nsample)

con = file(sprintf("%s.VarByCol.hap", prefix), "r", blocking = FALSE)
nsample = 0
while(length(line <- readLines(con, n = 1)) == 1) {
	nsample = nsample + 1
	line = strsplit(line, "", T)[[1]]
	line = as.integer(line)
	H[, nsample] = line
	if (nsample %% 100 == 0) cat(".")
}
cat("\n")
close(con)

cat("# haplotypes: ", ncol(H), "\n")
cat("# markers:    ", nrow(H), "\n")


# make genotypes
cat("Making genotypes\n")
G = matrix(integer(), nrow = nrow(H), ncol = ncol(H) / 2)
j = 1
for (i in 1:ncol(G)) {
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
	if (i %% 100 == 0) cat(".")
}
cat("\n")


# get positions
POS = as.numeric(readLines(sprintf("%s.pos", prefix), n = -1))


# get frequencies
cat("Getting frequencies\n")
MAC = rowSums(H)
MAF = MAC / ncol(H)
FRQ = MAF
i = which(MAC > (ncol(H) / 2))
if (length(i) > 0) {
	MAC[i] = ncol(H) - MAC[i]
	MAF[i] = 1 - MAF[i]
}


# function to parse identified sharers
parse.sharers = function(x) {
	if (length(x) == 1) {
		as.integer(strsplit(x, "|", T)[[1]])
	} else {
		lapply(strsplit(x, "|", T), as.integer)
	}
}


# function to address an individual's haplotypes in matrix
convert.gen2hap.sharers = function(i) {
	sort(c(i * 2 - 1, i * 2))
}


# make fk index
cat("Compiling rare variant index\n")
FKI = which(MAC > 1 & MAC <= fk.max)
cat("# rare variants: ", length(FKI), "\n")
FKI = lapply(FKI, function(i) {
	g = which(G[i, ] == 1)
	h = if (FRQ[i] <= 0.5) which(H[i, ] == 1) else which(H[i, ] == 0)
	h = intersect(h, convert.gen2hap.sharers(g))
	n = length(g)
	if (n < 2) return(NULL)
	data.table(index = i, 
						 pos = POS[i], 
						 frq = FRQ[i],
						 maf = MAF[i],
						 mac = MAC[i],
						 n.sharer = n, 
						 g.sharer = paste(as.character(g), collapse = "|"),
						 h.sharer = paste(as.character(h), collapse = "|"))
})
FKI = rbindlist(FKI)



# storing on disk
cat("Storing data on disk\n")

if (file.exists(sprintf("_bigmatrix.%s.H.bin",  prefix)) ||
		file.exists(sprintf("_bigmatrix.%s.H.desc",  prefix))) {
	unlink(c(sprintf("_bigmatrix.%s.H.bin",  prefix), 
					 sprintf("_bigmatrix.%s.H.desc",  prefix)))
}

tmp = big.matrix(nrow(H), ncol(H), type = "char", 
								 backingpath = getwd(),
								 backingfile =    sprintf("_bigmatrix.%s.H.bin",  prefix), 
								 descriptorfile = sprintf("_bigmatrix.%s.H.desc", prefix))

for (i in 1:ncol(H)) {
	tmp[, i] = H[, i]
	if (i %% 100 == 0) cat(".")
}
cat("\n")

if (file.exists(sprintf("_bigmatrix.%s.G.bin",  prefix)) ||
		file.exists(sprintf("_bigmatrix.%s.G.desc",  prefix))) {
	unlink(c(sprintf("_bigmatrix.%s.G.bin",  prefix), 
					 sprintf("_bigmatrix.%s.G.desc",  prefix)))
}

tmp = big.matrix(nrow(G), ncol(G), type = "char", 
								 backingpath = getwd(),
								 backingfile =    sprintf("_bigmatrix.%s.G.bin",  prefix), 
								 descriptorfile = sprintf("_bigmatrix.%s.G.desc", prefix))

for (i in 1:ncol(G)) {
	tmp[, i] = G[, i]
	if (i %% 100 == 0) cat(".")
}
cat("\n")

H = sprintf("%s.H", prefix)
G = sprintf("%s.G", prefix)



# function to load bigatrix data
load.bigmatrix = function(input.file) {
	attach.big.matrix(dget( sprintf("_bigmatrix.%s.desc", input.file)))
}


# save data
cat("Saving data\n")
save(H, G, 
		 POS, MAC, MAF, FKI, 
		 parse.sharers,
		 convert.gen2hap.sharers,
		 load.bigmatrix,
		 file = sprintf("%s.RData", prefix))


#
# write sharing pair table for IBD detection
#

cat("Distributing sharing pair tables ... \n")

rng = 1:nrow(FKI)
rng = split(rng, cut(rng, shapairs))

num.pairs = 0

for (k in 1:length(rng)) {
	cat(sprintf("%d of %d\n", k, length(rng)))
	
	ibd.name = sprintf("truth.%s.%04d.txt", prefix, k)
	
	ibd.file = file(ibd.name, open = "w")
	
	for (i in rng[[k]]) {
		idx = FKI$index[i]
		pos = FKI$pos[i]
		fk  = FKI$n.sharer[i]
		set = parse.sharers(FKI$h.sharer[i])
		set = combn(set, 2, simplify = F)
		for (pair in set) {
			cat(sprintf("%d %.8f %d %d %d\n", idx, pos, fk, pair[1], pair[2]), file = ibd.file)
			num.pairs = num.pairs + 1
		}
	}
	
	close(ibd.file)
}

cat("Number of pairs: ", num.pairs, "\n")








