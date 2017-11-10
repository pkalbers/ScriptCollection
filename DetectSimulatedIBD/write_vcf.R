#
# write VCF file for Beagle IBD
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]

load(sprintf("%s.RData", prefix))


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


cat("Writing VCF files\n")

H = load.bigmatrix(sprintf("%s.H", prefix))
P = load.bigmatrix(sprintf("%s.P", prefix))

con.h = file(sprintf("%s.H.vcf", prefix), open = "w")
con.p = file(sprintf("%s.P.vcf", prefix), open = "w")

head = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
indv = sprintf("ID%05d", 1:(ncol(H) / 2))
head = paste(head, paste(indv, collapse = "\t"), sep = "\t")

cat(head, "\n", sep = "", file = con.h)
cat(head, "\n", sep = "", file = con.p)

for (i in 1:nrow(H)) {
	h = H[i, ]
	p = P[i, ]
	
	a0 = sample(c("A", "G", "C", "T"), 1)
	a1 = sample(setdiff(c("A", "G", "C", "T"), a0), 1)
	
	inf = sprintf("1\t%d\tsnp_id_%08d\t%s\t%s\t100\tPASS\t.\tGT\t", pos[i], pos[i], a0, a1)
	
	hap = matrix(h, nrow = ncol(H) / 2, ncol = 2, byrow = T)
	hap = apply(hap, 1, function(x) sprintf("%d|%d", x[1], x[2]) )
	hap = paste(hap, collapse = "\t")
	
	phap = matrix(p, nrow = ncol(H) / 2, ncol = 2, byrow = T)
	phap = apply(phap, 1, function(x) sprintf("%d|%d", x[1], x[2]) )
	phap = paste(phap, collapse = "\t")
	
	cat(inf, hap,  "\n", sep = "", file = con.h)
	cat(inf, phap, "\n", sep = "", file = con.p)
	if (i %% 10000 == 0) cat(".")
}
cat("\n")

close(con.h)
close(con.p)

cat("DONE\n")




