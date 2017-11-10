#
# load simulated data into R
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

fk.max = as.numeric(args[1]) # rare variant max count
prefix = args[2]             # file prefix
shapairs = args[3]           # number of distributed sharing pairs files



# get simulated haplotypes
cat("Reading haplotypes\n")
tmp = readLines(sprintf("%s.VarByCol.hap", prefix), n = -1)
tmp = strsplit(tmp, "", T)
tmp = lapply(tmp, as.integer)
H = matrix(0, nrow = length(tmp[[1]]), ncol = length(tmp))
for (i in 1:length(tmp)) {
	H[, i] = tmp[[i]]
}
cat("# haplotypes: ", ncol(H), "\n")
cat("# markers:    ", nrow(H), "\n")


# make genotypes
cat("Making genotypes\n")
G = matrix(0, nrow = nrow(H), ncol = ncol(H) / 2)
j = 1
for (i in 1:ncol(G)) {
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
}


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
	}
	lapply(strsplit(x, "|", T), as.integer)
}


# function to address an individual's haplotypes in matrix
convert.gen2hap.sharers = function(i) {
	sort(c(i * 2 - 1, i * 2))
}


# make fk index
cat("Compiling rare variant index\n")
FKI = which(MAC > 1 & MAC <= fk.max)
FKI = lapply(FKI, function(i) {
	g = which(G[i, ] == 1)
	h = if (frq[i] <= 0.5) which(H[i, ] == 1) else which(H[i, ] == 0)
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


H = as.big.matrix(H, type="integer",
									backingpath = getwd(),
									backingfile =  sprintf("_bigmatrix.%s.bin", prefix), 
									descriptorfile = sprintf("_bigmatrix.%s.desc", prefix))

G = as.big.matrix(G, type="integer",
									backingpath = getwd(),
									backingfile =  sprintf("_bigmatrix.%s.bin", prefix), 
									descriptorfile = sprintf("_bigmatrix.%s.desc", prefix))


# function to load bigatrix data
load.bigmatrix = function(input.file) {
	attach.big.matrix(dget( sprintf("bigmatrix.%s.desc", input.file)))
}

# save data
cat("Saving data\n")
save(H, G, 
		 POS, MAC, MAF, FKI, 
		 parse.sharers,
		 get.indv.haps,
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
	
	ibd.name = sprintf("sharingpairs.%s.%04d.txt", prefix, k)
	
	ibd.file = file(ibd.name, open = "w")
	
	for (i in rng[[k]]) {
		idx = FKI$index[i]
		pos = FKI$pos[i]
		set = parse.sharers(FKI$h.sharer[i])
		set = combn(set, 2, simplify = F)
		for (pair in set) {
			cat(sprintf("%d %.8f %d %d\n", idx, pos, pair[1], pair[2]), file = ibd.file)
			num.pairs = num.pairs + 1
		}
	}
	
	close(ibd.file)
}

cat("Number of pairs: ", num.pairs, "\n")








