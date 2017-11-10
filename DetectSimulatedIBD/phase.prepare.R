#
# write GEN + sample files for shapeit
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]

load(sprintf("%s.RData", prefix))


p = round(POS)

i = anyDuplicated(p)
while (i != 0) {
	if (i == 1) {
		p[i] = p[i] - 1
	} else {
		p[i] = p[i] + 1
	}
	i = anyDuplicated(p)
}


cat("Writing GEN file\n")

con = file(sprintf("%s.gen", prefix), open = "w")

G = load.bigmatrix(G)

for (i in 1:nrow(G)) {
	g = G[i, ]
	varinf = sprintf("0 SNP%08d %d A B ", i, p[i])
	genvec = rep("", ncol(G))
	for (j in 1:ncol(G)) {
		if (g[j] == 0) { genvec[j] = "1 0 0"; next }
		if (g[j] == 1) { genvec[j] = "0 1 0"; next }
		if (g[j] == 2) { genvec[j] = "0 0 1"; next }
	}
	genvec = paste(genvec, collapse = " ")
	cat(varinf, genvec, "\n", sep = "", file = con)
	if (i %% 10000 == 0) cat(".")
}
cat("\n")

close(con)


cat("Writing Sample file\n")
con = file(sprintf("%s.sample", prefix), open = "w")
cat("ID_1 ID_2 missing\n", file = con)
cat("0 0 0\n", file = con)
for (i in 1:ncol(G)) {
	cat(sprintf("IND%04d IND%04d 0\n", i, i), file = con)
	if (i %% 100 == 0) cat(".")
}
cat("\n")
close(con)





