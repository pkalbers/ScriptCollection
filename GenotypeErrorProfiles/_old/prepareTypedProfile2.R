#
# GT error profile, assuming that non-contained sites are homref, after matching with freq from 1000G
#

args = commandArgs(T)

vcf.file = args[1]
frq.file = args[2]



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



typed = data.frame(key = vcf$KEY, 
									 ref = vcf$REF, 
									 alt = vcf$ALT, 
									 frq.ref = NA,
									 frq.alt = NA,
									 gt = vcf$GT, 
									 stringsAsFactors = F)

rownames(typed) = vcf$KEY




cat("Compiling genotype frequencies ... ")

ref.frq = rep(NA, nrow(vcf))
alt.frq = rep(NA, nrow(vcf))
names(ref.frq) = vcf$KEY
names(alt.frq) = vcf$KEY

sub.vcf = vcf[key, ]
sub.frq = frq[key, ]

inv = which(sub.vcf$REF == sub.frq$alt & sub.vcf$ALT == sub.frq$ref)
if (length(inv) > 0) {
	cat(sprintf("(inverted alleles: %d) ", length(inv)))
	tmp = sub.frq$ref.frq[inv]
	sub.frq$ref.frq[inv] = sub.frq$alt.frq[inv]
	sub.frq$alt.frq[inv] = tmp
}

a = which(sub.vcf$REF != sub.frq$ref)
b = which(sub.vcf$ALT != sub.frq$alt)
del = sort(unique(c(a, b)))
if (length(del) > 0) {
	cat(sprintf("(non-matching alleles: %d) ", length(del)))
	sub.frq = sub.frq[-del, ]
}

key = rownames(sub.frq)

ref.frq[ key ] = sub.frq$ref.frq
alt.frq[ key ] = sub.frq$alt.frq

typed$frq.ref = ref.frq
typed$frq.alt = alt.frq

cat("OK\n")





cat("Saving ... ")
save(typed, file = "typed.RData")
cat("DONE\n")








