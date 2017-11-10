#
# Load real 1KG
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)


fk.max = 25 # as.numeric(args[1]) # rare variant max count
shapairs = 1024 # as.numeric(args[3]) # number of distributed sharing pairs files

prefix = "1000GP_Phase3_chr20"



cat("Reading haplotypes\n")

con = file(sprintf("%s.hap", prefix), "r", blocking = FALSE)
nmarker = 0
nsample = 0
while(length(line <- readLines(con, n = 1)) == 1) {
	x = nchar(line)
	if (x <= 1) {
		break
	}
	if (nsample == 0) {
		nsample = x
	} else {
		if (nsample != nchar(line)) stop("Different number of markers")
	}
	nmarker = nmarker + 1
	if (nmarker %% 100000 == 0) cat(".")
}
cat("\n")
close(con)


nsample = ceiling(nsample / 2)


H = matrix(integer(), nrow = nmarker, ncol = nsample)

con = file(sprintf("%s.hap", prefix), "r", blocking = FALSE)
nmarker = 0
while(length(line <- readLines(con, n = 1)) == 1) {
	x = nchar(line)
	if (x <= 1) {
		break
	}
	nmarker = nmarker + 1
	line = strsplit(line, " ", T)[[1]]
	line = as.integer(line)
	H[nmarker, ] = line
	if (nmarker %% 100000 == 0) cat(".")
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
POS = read.table(sprintf("%s.legend", prefix), header = T, stringsAsFactors = F)
POS = POS$position


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




# save data
cat("Saving data\n")

save(H, file = sprintf("%s.H.RData", prefix))
save(G, file = sprintf("%s.G.RData", prefix))

save(POS, MAC, MAF, FKI, 
		 parse.sharers,
		 convert.gen2hap.sharers,
		 file = sprintf("%s.RData", prefix))





dir.create("pairs", showWarnings = F)
setwd("pairs")


cat("Distributing sharing pair tables ... \n")

rng = 1:nrow(FKI)
rng = split(rng, cut(rng, shapairs))

num.pairs = 0

for (k in 1:length(rng)) {
	cat(sprintf("%d of %d  ", k, length(rng)))
	
	pair.name = sprintf("pairs.%s.%04d.RData", prefix, k)
	
	pair = NULL
	
	for (i in rng[[k]]) {
		
		g = parse.sharers(FKI$g.sharer[i])
		h = parse.sharers(FKI$h.sharer[i])
		
		g = combn(g, 2, simplify = F)
		h = combn(h, 2, simplify = F)
		
		tmp = data.table(index = FKI$index[i],
										 position = FKI$pos[i],
										 fk = FKI$n.sharer[i],
										 g0 = sapply(g, function(x) x[1]),
										 g1 = sapply(g, function(x) x[2]),
										 h0 = sapply(h, function(x) x[1]),
										 h1 = sapply(h, function(x) x[2]))
		
		pair = rbind(pair, tmp)
		
		num.pairs = num.pairs + nrow(tmp)
	}
	
	save(pair, file = pair.name)
	cat("\n")
}

cat("Number of pairs: ", num.pairs, "\n")







