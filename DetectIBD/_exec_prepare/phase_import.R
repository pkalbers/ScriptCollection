#
# read phased haplotypes
#

library(data.table)

#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20.GenErr_1000G" #args[1]
phap.file = "phased.OutOfAfricaHapMap20.GenErr_1000G.haps" # args[2] # phased file


load(sprintf("data.%s.RData", prefix))
load(sprintf("data.%s.H.RData", prefix)) # haplotypes





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


save(P, file = sprintf("data.%s.P.RData", prefix))


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

a = FKI$g0 * 2 - 1
b = FKI$g0 * 2
coord.a = as.matrix(data.frame(x = FKI$index, y = a))
coord.b = as.matrix(data.frame(x = FKI$index, y = b))
coord.a = P[ coord.a ]
coord.b = P[ coord.b ]

if (any(coord.a == coord.b)) stop("???")
FKI$p0 = (a * coord.a) + (b * coord.b)

a = FKI$g1 * 2 - 1
b = FKI$g1 * 2
coord.a = as.matrix(data.frame(x = FKI$index, y = a))
coord.b = as.matrix(data.frame(x = FKI$index, y = b))
coord.a = P[ coord.a ]
coord.b = P[ coord.b ]
if (any(coord.a == coord.b)) stop("???")
FKI$p1 = (a * coord.a) + (b * coord.b)



# save data
cat("Saving data\n")

if ("applied.error.table" %in% ls() && "applied.error.matrix" %in% ls()) {
	
	save(POS, CPOS, BRK, CBRK, 
			 AAC, AAF, MAC, MAF, 
			 GC0, GC1, GC2, GF0, GF1, GF2,
			 FKI,
			 applied.error.table, applied.error.matrix,
			 file = sprintf("data.%s.RData", prefix))
	
} else {
	
	save(POS, CPOS, BRK, CBRK, 
			 AAC, AAF, MAC, MAF, 
			 GC0, GC1, GC2, GF0, GF1, GF2,
			 FKI,
			 file = sprintf("data.%s.RData", prefix))
	
}

cat("OK\n")




