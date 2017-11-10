#
# prepare GT truth data
#

args = commandArgs(T)

vcf.file = args[1]
bed.file = args[2]
frq.file = args[3]
out.file = args[4]



cat("Loading data ... ")
vcf = read.table(vcf.file, header = F, stringsAsFactors = F)
bed = read.table(bed.file, header = F, stringsAsFactors = F)
load(frq.file)
cat("OK\n")



header = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
names(vcf) = header

header = c("CHROM", "BEG", "END")
names(bed) = header



cat("Adjusting chromosomes ... ")

vcf$CHROM = as.character(vcf$CHROM)
if (any(grepl("^chr.+$", unique(vcf$CHROM)))) {
	vcf$CHROM = sub("^chr(.+)$", "\\1", vcf$CHROM)
}
del = grep("^[0-9]+$", vcf$CHROM, invert = T)
if (length(del) > 0) {
	vcf = vcf[-del, ]
}

bed$CHROM = as.character(bed$CHROM)
if (any(grepl("^chr.+$", unique(bed$CHROM)))) {
	bed$CHROM = sub("^chr(.+)$", "\\1", bed$CHROM)
}
del = grep("^[0-9]+$", bed$CHROM, invert = T)
if (length(del) > 0) {
	bed = bed[-del, ]
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



cat("Compiling true genotypes ... ")

sub.vcf = vcf[key, ]
sub.frq = frq[key, ]

a = which(sub.vcf$REF != sub.frq$ref)
b = which(sub.vcf$ALT != sub.frq$alt)
del = sort(unique(c(a, b)))
if (length(del) > 0) {
	cat(sprintf("(removing non-matching alleles: %d) ", length(del)))
	sub.vcf = sub.vcf[-del, ]
}

truth = sub.vcf$GT
names(truth) = sub.vcf$KEY

cat("OK\n")




cat("Compiling masked true genotypes ")

sub.vcf = split(sub.vcf, sub.vcf$CHROM)
sub.bed = split(bed, bed$CHROM)

tags = intersect(names(sub.bed), names(sub.vcf))

truth.mask = c()

for (tag in tags) {
	cat(".")
	cv = sub.vcf[[tag]]
	cb = sub.bed[[tag]]
	
	end = max(max(cb$END), max(cv$POS))
	
	x = rep(FALSE, end)
	for (i in 1:nrow(cb)) {
		x[(cb$BEG[i]):(cb$END[i])] = TRUE
	}
	y = which(x[cv$POS])
	
	if (length(y) > 0) {
		truth.mask = c(truth.mask, cv$KEY[y])	
	}
}

cat("OK\n")



cat("Assuming confident hom-ref genotypes ")

sub.frq = frq[-(which(frq$key %in% key)), ]
sub.frq = split(sub.frq, sub.frq$chr)
#sub.bed = split(bed, bed$CHROM)

tags = intersect(names(sub.bed), names(sub.frq))

conf.hom.ref = c()

for (tag in tags) {
	cat(".")
	cf = sub.frq[[tag]]
	cb = sub.bed[[tag]]
	
	end = max(max(cb$END), max(cf$pos))
	
	x = rep(FALSE, end)
	for (i in 1:nrow(cb)) {
		x[(cb$BEG[i]):(cb$END[i])] = TRUE
	}
	y = which(x[cf$pos])

	if (length(y) > 0) {
		conf.hom.ref = c(conf.hom.ref, cf$key[y])	
	}
}

cat(" OK\n")



cat("Saving ... ")
save(truth, truth.mask, conf.hom.ref, file = out.file)
cat("DONE\n")




