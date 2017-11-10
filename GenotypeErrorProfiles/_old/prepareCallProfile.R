#
# GT error profile, assuming that non-contained sites are homref, after matching with freq from 1000G
#

args = commandArgs(T)

vcf.file = args[1]
frq.file = args[2]
out.file = args[3]



cat("Loading data ... ")
vcf = read.table(vcf.file, header = F, stringsAsFactors = F)
load(frq.file)
cat("OK\n")



header = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
names(vcf) = header



cat("Adjusting chromosomes ... ")
vcf$CHROM = as.character(vcf$CHROM)
if (any(grepl("^chr.+$", unique(vcf$CHROM)))) {
	vcf$CHROM = sub("^chr(.+)$", "\\1", vcf$CHROM)
}
del = grep("^[0-9]+$", vcf$CHROM, invert = T)
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
cat("OK\n")



cat("Retaining SNPs only ... ")
a = which(nchar(vcf$REF) > 1)
b = which(nchar(vcf$ALT) > 1)
c = which(vcf$REF == ".")
d = which(vcf$ALT == ".")
del = sort(unique(c(a,b,c,d)))
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
cat("OK\n")



cat("Adjusting genotypes ... ")
del = unique(c(grep("^\\.[\\/\\|].+$", vcf$SAMPLE), grep("^[01\\.][\\/\\|]\\..*$", vcf$SAMPLE)))
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
del = grep("^[01][\\/\\|][01].*", vcf$SAMPLE, invert = T)
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
vcf$GT = as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\1", vcf$SAMPLE)) + as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\2", vcf$SAMPLE))
cat("OK\n")



cat("Generating keys ... ")
vcf$KEY = paste(vcf$CHROM, as.character(vcf$POS), sep = "_")
del = which(duplicated(vcf$KEY))
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
rownames(vcf) = vcf$KEY
cat("OK\n")



key = intersect(vcf$KEY, frq$key)
cat("Variants matching to known frequencies:", length(key), sprintf("(%.3f%%)", length(key) / nrow(vcf) * 100), "\n")


cat("Compiling profile ... ")

vcf = vcf[key, ]
sub.frq = frq[key, ]


p = cbind(sub.frq, 
					call.gt = vcf$GT,
					call.is.typed = T,
					true.gt = NA,
					true.is.typed = NA,
					true.is.confident = NA)


a = which(p$ref != vcf$REF)
b = which(p$alt != vcf$ALT)
del = sort(unique(c(a, b)))
if (length(del) > 0) {
	cat(sprintf("(removing non-matching alleles: %d) ", length(del)))
	p = p[-del, ]
}


if (nrow(frq) != length(key)) {
	sub.frq = frq[-(which(frq$key %in% key)), ]
	
	q = cbind(sub.frq, 
						call.gt = 0,
						call.is.typed = F,
						true.gt = NA,
						true.is.typed = NA,
						true.is.confident = NA)
	
	
	p = rbind(p, q)
}

p$key = NULL

cat("OK\n")


cat("Saving ... ")
save(p, file = out.file)
cat("DONE\n")








