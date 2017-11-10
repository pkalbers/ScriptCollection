#
# read phased haplotypes
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]             # file prefix
phap.file = args[2]          # phased file


load(sprintf("%s.RData", prefix))


H = load.bigmatrix(sprintf("%s.H", prefix))


# get phased haplotypes
cat("Reading phased haplotypes ... \n")

P = matrix(integer(), nrow = nrow(H), ncol = ncol(H))

con = file(phap.file, "r", blocking = FALSE)
i = 0
n = ncol(H)
while(length(line <- readLines(con, n = 1)) == 1) {
	i = i + 1
	line = strsplit(line, " ", T)[[1]]
	line = line[6:length(line)]
	if (length(line) != n) stop("Different sample size on line: ", i)
	line = as.integer(line)
	P[i, ] = line
	if (i %% 10000 == 0) cat(i, "of", nrow(H), "\n")
}
cat("Done\n")
close(con)


# consider one full switch
# cat("Switching haplotypes\n")
# j = 1
# for (i in 1:(ncol(H) / 2)) {
# 	h0 = H[, j]
# 	h1 = H[, j+1]
# 	
# 	p0 = P[, j]
# 	p1 = P[, j+1]
# 	
# 	a = sum(abs(h0 - p0)) + sum(abs(h1 - p1))
# 	b = sum(abs(h0 - p1)) + sum(abs(h1 - p0))
# 	
# 	if (b < a) {
# 		P[, j]   = p1
# 		P[, j+1] = p0
# 	}
# 	
# 	j = j + 2
# 	if (i %% 100 == 0) cat(".")
# }
# cat("\n")



# check frequencies
cat("Checking frequencies\n")
for (i in 1:nrow(H)) {
	if (sum(H[i, ]) != sum(P[i, ])) {
		stop("Unequal frequencies!")
	}
	if (i %% 10000 == 0) cat(i, "of", nrow(H), "\n")
}
cat("Done\n")



# completing fk index
cat("Adding to rare variant index\n")
FKI$p.sharer = ""
for (i in 1:nrow(FKI)) {
	idx = FKI$index[i]
	set = FKI$g.sharer[i]
	set = parse.sharers(set)
	set = convert.gen2hap.sharers(set)
	det = if (FKI$frq[i] <= 0.5) which(P[idx, set] == 1) else which(P[idx, set] == 0)
	p = set[det]
	FKI$p.sharer[i] = paste(as.character(p), collapse = "|")
	if (i %% 5000 == 0) cat(i, "of", nrow(FKI), "\n")
}



# storing on disk
cat("Storing data on disk\n")

if (file.exists(sprintf("_bigmatrix.%s.P.bin",  prefix)) ||
		file.exists(sprintf("_bigmatrix.%s.P.desc",  prefix))) {
	unlink(c(sprintf("_bigmatrix.%s.P.bin",  prefix), 
					 sprintf("_bigmatrix.%s.P.desc",  prefix)))
}

tmp = big.matrix(nrow(P), ncol(P), type = "char", 
								 backingpath = getwd(),
								 backingfile =    sprintf("_bigmatrix.%s.P.bin",  prefix), 
								 descriptorfile = sprintf("_bigmatrix.%s.P.desc", prefix))

for (i in 1:ncol(P)) {
	tmp[, i] = P[, i]
	if (i %% 100 == 0) cat(".")
}
cat("\n")


H = sprintf("%s.H", prefix)
G = sprintf("%s.G", prefix)
P = sprintf("%s.P", prefix)


# save data
cat("Saving data\n")
save(H, G, P,
		 POS, MAC, MAF, FKI, 
		 parse.sharers,
		 convert.gen2hap.sharers,
		 load.bigmatrix,
		 applied.err.mat, applied.err.frq, err.data, err.func, ### <--
		 file = sprintf("%s.RData", prefix))




