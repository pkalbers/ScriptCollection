
library(bigmemory)
library(data.table)
library(digest)
options(bigmemory.typecast.warning=FALSE)


args = commandArgs(TRUE)

hap.file = sprintf("%s.hap", args[1])
pos.file = sprintf("%s.pos", args[1])
gen.file = sprintf("%s.gen", args[1])
shh.file = sprintf("%s.shh", args[1])

f.level = 25
chrom.length = 50e06


# haplotypes / genotypes, individuals

lines = readLines(hap.file, n = -1)

H = lapply(strsplit(lines, "", TRUE), as.integer)
H = as.data.frame(H)
names(H) = as.character(1:length(H))

G  = matrix(NA, length(H[[1]]), length(H) / 2)

I = c()

i = 0
k = 0
while (i < length(H)) {
	k = k + 1
	
	i = i + 1
	a = as.character(i)
	i = i + 1
	b = as.character(i)

	G[, k]  = H[[a]] + H[[b]]
	
	I = c(I, sprintf("%s %s", a, b))
}

colnames(G) = I

names(I) = sprintf("INDV%05d", 1:length(I))
J = names(I)
names(J) = I

K = rainbow(length(I))
K = K[order(sapply(I, digest))]
names(K) = I


H = as.matrix(H)
H = as.big.matrix(H, type="char",
									backingpath = getwd(),
									backingfile =  sprintf("bigmatrix.%s.bin", hap.file), 
									descriptorfile = sprintf("bigmatrix.%s.desc", hap.file))

G = as.big.matrix(G, type="char",
									backingpath = getwd(),
									backingfile =  sprintf("bigmatrix.%s.bin", gen.file), 
									descriptorfile = sprintf("bigmatrix.%s.desc", gen.file))


# positions

lines = readLines(pos.file, n = -1)
pos = strsplit(lines, " ", TRUE)[[1]]
pos = as.numeric(pos)
pos = round(pos * chrom.length)

dup = which(duplicated(pos))
while (length(dup) != 0) {
	pos[dup] = pos[dup] + 1
	dup = which(duplicated(pos))
	cat(".\n")
}
if (!identical(pos, sort(pos))) stop("Positions!")

M = data.frame(position = pos, af0 = NA, af1 = NA, ac0 = NA, ac1 = NA, f = 0, stringsAsFactors = FALSE)

M$ac0 = apply(H[, ], 1, function(x) length(which(x == 0)))
M$ac1 = apply(H[, ], 1, function(x) length(which(x == 1)))
M$af0 = M$ac0 / ncol(H)
M$af1 = M$ac1 / ncol(H)


# sharing

S = matrix(0, nrow(G), ncol(G), dimnames = list(NULL, as.vector(I)))

for (i in 1:nrow(M)) {
	f = min(M$ac0[i], M$ac1[i])
	if (f > 1 && f <= f.level) {
		col = which(G[i, ] == 1)
		lev = length(col)
		if (lev > 1 && lev <= f.level) {
			M$f[i] = lev
			S[i, col] = 1
		}
	}
}

S = as.big.matrix(S, type="char",
									backingpath = getwd(),
									backingfile =  sprintf("bigmatrix.%s.bin", shh.file), 
									descriptorfile = sprintf("bigmatrix.%s.desc", shh.file))




load.bigmatrix = function(input.file) {
	attach.big.matrix(dget( sprintf("bigmatrix.%s.desc", input.file)))
}

save(H, G, I, J, K, M, S, 
		 hap.file, pos.file, gen.file, shh.file, 
		 load.bigmatrix, 
		 file = sprintf("data.%s.RData", args[1]))

