#
# prepare GT truth data
#

library(data.table)

args = commandArgs(T)

vcf.file = args[1]
bed.file = args[2]
frq.file = args[3]
prefix = args[4]



cat("Loading data ... ")
vcf = read.table(vcf.file, header = F, stringsAsFactors = F)
bed = read.table(bed.file, header = F, stringsAsFactors = F)
load(frq.file)
cat("OK\n")


cat(sprintf("# Variants: %d\n", nrow(vcf)))


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


cat(sprintf("# Variants: %d\n", nrow(vcf)))


cat("Retaining SNPs only ... ")
a = which(nchar(vcf$REF) > 1)
b = which(nchar(vcf$ALT) > 1)
del = sort(unique(c(a,b)))
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
cat("OK\n")


cat(sprintf("# Variants: %d\n", nrow(vcf)))


cat("Parsing genotypes ... ")
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


cat(sprintf("# Variants: %d\n", nrow(vcf)))


cat("Generating variant IDs ... ")
vcf$KEY = paste(vcf$CHROM, as.character(vcf$POS), sep = "_")
del = which(duplicated(vcf$KEY))
if (length(del) > 0) {
	vcf = vcf[-del, ]
}
rownames(vcf) = vcf$KEY
cat("OK\n")


cat(sprintf("# Variants: %d\n", nrow(vcf)))


# splitting by chromosome
cat("\n")

cat("Preparing data by chromosome ... ")

vcf = split(vcf, vcf$CHROM)
bed = split(bed, bed$CHROM)
frq = split(frq, frq$chr)

chr = as.character(sort(as.integer(intersect(names(vcf), intersect(names(bed), names(frq))))))

cat("OK\n")



data = list()

for (tag in chr) {
	cat("# Chromosome:", tag, "\n")
	
	v = vcf[[tag]]
	f = frq[[tag]]
	b = bed[[tag]]
	
	rng = rep(FALSE, max(c(max(b$END), max(v$POS), max(f$pos))))
	
	for (i in 1:nrow(b)) {
		seg = (b$BEG[i]) : (b$END[i])
		rng[seg] = TRUE
	}
	
	vcf.confident = rng[ v$POS ]
	frq.confident = rng[ f$pos ]
	
	keep = which(vcf.confident)
	cat(sprintf("  Confident variants in VCF file: %d of %d (%.1f%%)\n", length(keep), nrow(v), length(keep) / nrow(v) * 100))
	v = v[keep, ]
	
	keep = which(frq.confident)
	cat(sprintf("  Confident variants in FRQ file: %d of %d (%.1f%%)\n", length(keep), nrow(f), length(keep) / nrow(f) * 100))
	f = f[keep, ]
	
	
	tmp = data.frame(chr = as.integer(tag), 
									 pos = v$POS, 
									 ref = v$REF, 
									 alt = v$ALT,
									 frq.ref = NA,
									 frq.alt = NA,
									 true.gt  = v$GT, 
									 true.is.assumed = F,
									 stringsAsFactors = F)
	
	rownames(tmp) = v$KEY
	
	key = intersect(v$KEY, f$key)
	cat("  Variants matching to known frequencies:", length(key), sprintf("(%.1f%%)", length(key) / nrow(v) * 100), "\n")
	
	tmp[key, ]$frq.ref = f[key, ]$ref.frq
	tmp[key, ]$frq.alt = f[key, ]$alt.frq
	
	del = which(tmp[key, ]$ref != f[key, ]$ref | tmp[key, ]$alt != f[key, ]$alt)
	cat("  Mismatched alleles at known frequencies:", length(del), "\n")
	if (length(del) > 0) {
		del = which(v$KEY %in% key[del])
		tmp = tmp[-del, ]
	}
	
	rownames(tmp) = NULL
	tru = tmp
	
	key = setdiff(f$key, v$KEY)
	cat("  Assuming '0' genotypes at known sites:", length(key), sprintf("(%.1f%%)", length(key) / nrow(f) * 100), "\n")
	
	f = f[key, ]
	
	tmp = data.frame(chr = as.integer(tag), 
									 pos = f$pos, 
									 ref = f$ref, 
									 alt = f$alt,
									 frq.ref = f$ref.frq,
									 frq.alt = f$alt.frq,
									 true.gt  = 0, 
									 true.is.assumed = T,
									 stringsAsFactors = F)
	
	tru = rbind(tru, tmp)
	
	cat(sprintf("  # Variants retained: %d\n", nrow(data[[tag]])))
	
	data[[tag]] = tru
}


cat(sprintf("# Overall variants retained: %d\n", sum(sapply(data, nrow))))

conf = bed

cat("Saving ... ")
save(data, conf, file = sprintf("truth.%s.RData", prefix))
cat("DONE\n")




####

load("truth.NA12878.RData")


sum(sapply(data, function(d) {
	length(which(!is.na(d$frq.ref) & !is.na(d$frq.alt) & !d$true.is.assumed))
}))

sum(sapply(data, function(d) {
	length(which(!d$true.is.assumed))
}))







