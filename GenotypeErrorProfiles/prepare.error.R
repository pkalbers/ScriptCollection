#
# GT error profile
#

args = commandArgs(T)

vcf.file = args[1]
tru.file = args[2]
prefix = args[3]


cat("Loading data ... ")
vcf = read.table(vcf.file, header = F, stringsAsFactors = F)
load(tru.file)
cat("OK\n")


cat(sprintf("# Variants: %d\n", nrow(vcf)))


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


cat("\n")


vcf = split(vcf, vcf$CHROM)

chr = as.character(sort(as.integer(intersect(names(vcf), intersect(names(conf), names(data))))))

p = NULL

for (tag in chr) {
	cat("# Chromosome:", tag, "\n")
	
	v = vcf[[tag]]
	b = conf[[tag]]
	d = data[[tag]]
	
	vcf.key = paste(v$CHROM, as.character(v$POS), sep = "_")
	dat.key = paste(as.character(d$chr), as.character(d$pos), sep = "_")
	
	if (anyDuplicated(vcf.key)) {
		del = which(duplicated(vcf.key))
		vcf.key = vcf.key[-del]
		v = v[-del, ]
	}
	
	if (anyDuplicated(dat.key)) {
		del = which(duplicated(dat.key))
		dat.key = dat.key[-del]
		d = d[-del, ]
	}
	
	key = intersect(vcf.key, dat.key)
	cat("  Variants matching to truth positions:", length(key), sprintf("(%.1f%%)", length(key) / nrow(v) * 100), "\n")
	
	d$call.gt = 0
	d$call.is.assumed = TRUE
	
	if (length(key) > 0) {
		
		vcf.idx = match(key, vcf.key)
		dat.idx = match(key, dat.key)
		
		d$call.gt[dat.idx] = v$GT[vcf.idx]
		d$call.is.assumed[dat.idx] = FALSE
		
		del0 = which(d$ref[dat.idx] != v$REF[vcf.idx])
		del1 = which(d$alt[dat.idx] != v$ALT[vcf.idx])
		del = union(del0, del1)
		if (length(del) > 0) {
			cat("  Removing mismatched alleles:", length(del), sprintf("(%.1f%%)", length(del) / length(key) * 100), "\n")
			del = dat.idx[del]
			dat.key = dat.key[-del]
			d = d[-del, ]
		}
		
	}
	
	key = setdiff(vcf.key, dat.key)
	cat("  Remaining variants in VCF file:", length(key), sprintf("(%.1f%%)", length(key) / nrow(v) * 100), "\n")
	
	if (length(key) > 1) {
		
		v = v[match(key, vcf.key), ]
		
		rng = rep(FALSE, max(c(max(b$END), max(v$POS))))
		
		for (i in 1:nrow(b)) {
			seg = (b$BEG[i]) : (b$END[i])
			rng[seg] = TRUE
		}
		
		vcf.confident = rng[ v$POS ]
		
		keep = which(vcf.confident)
		cat("  Matching to confident regions:", length(keep), sprintf("(%.1f%%)", length(keep) / nrow(v) * 100), "\n")
		
		if (length(keep) > 1) {
			
			v = v[keep, ]
			
			tmp = data.frame(chr = as.integer(tag), 
											 pos = v$POS, 
											 ref = v$REF, 
											 alt = v$ALT, 
											 frq.ref = NA, 
											 frq.alt = NA, 
											 true.gt = 0, 
											 true.is.assumed = T, 
											 call.gt = v$GT, 
											 call.is.assumed = F,
											 stringsAsFactors = F)
			
			d = rbind(d, tmp)
		}
	}
	
	
	del = which(d$call.is.assumed)
	if (length(del) > 0) {
		d = d[-del, ]
	}
	
	
	p = rbind(p, d)
}


cat("\n")

cat(sprintf("# Variants: %d\n", nrow(p)))


cat("Saving ... ")

save(p, file = sprintf("profile.%s.RData", prefix))

cat("DONE\n")








