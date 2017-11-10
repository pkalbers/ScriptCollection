#
# GT error profile
#

args = commandArgs(T)

ref.file = args[1] # "../source_illumina/hg19/8.0.1/NA12878/NA12878.vcf.gz"
alt.file = args[2] # "vcfBeta-NA12878-200-37-ASM.vcf.bz2" # "NA12878_Affymetrix_Axiom_DB_2010_v4_b37.recode.vcf"
bed.file = args[3] # "../source_illumina/hg19/8.0.1/NA12878/ConfidentRegions.bed"
frq.file = args[4] # "/Volumes/Sirius/GenotypeErrorProfiles/VCF/source_1000g/1000g/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.frq"
out.file = args[5]

cat(ref.file, "\n")
cat(alt.file, "\n")
cat(bed.file, "\n")
cat(frq.file, "\n")
cat(out.file, "\n\n")


cat("Loading data ... ")
ref = read.table(ref.file, header = F, stringsAsFactors = F)
alt = read.table(alt.file, header = F, stringsAsFactors = F)
bed = read.table(bed.file, header = F, stringsAsFactors = F)
frq = readLines(frq.file)
cat("OK\n")



cat("Parsing frequencies ... ")
frq = frq[-1]
frq = strsplit(frq, "\t")

tag = lapply(frq, function(x) {
	paste(x[1], x[2], sep = "_")
})
tag = unlist(tag)

names(frq) = tag

maf = sapply(frq, function(x) {
	x = x[-(1:4)]
	if (length(x) > 2) return(NULL)
	val = as.numeric(sub("^([A-Z]):(.+)$", "\\2", x))
	min(val)
})
maf = unlist(maf)
cat("OK\n")



cat("Adjusting headers ... ")
header = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
names(ref) = header
names(alt) = header

header = c("CHROM", "BEG", "END")
names(bed) = header
cat("OK\n")



cat("Adjusting chromosomes ... ")
ref$CHROM = as.character(ref$CHROM)
alt$CHROM = as.character(alt$CHROM)

if (any(grepl("^chr.+$", unique(ref$CHROM)))) {
	ref$CHROM = sub("^chr(.+)$", "\\1", ref$CHROM)
}

if (any(grepl("^chr.+$", unique(alt$CHROM)))) {
	alt$CHROM = sub("^chr(.+)$", "\\1", alt$CHROM)
}

del = grep("^[0-9]+$", ref$CHROM, invert = T)
if (length(del) > 0) {
	ref = ref[-del, ]
}

del = grep("^[0-9]+$", alt$CHROM, invert = T)
if (length(del) > 0) {
	alt = alt[-del, ]
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
a = which(nchar(ref$REF) > 1)
b = which(nchar(ref$ALT) > 1)
c = which(ref$REF == ".")
d = which(ref$ALT == ".")
del = sort(unique(c(a,b,c,d)))
if (length(del) > 0) {
	ref = ref[-del, ]
}

a = which(nchar(alt$REF) > 1)
b = which(nchar(alt$ALT) > 1)
c = which(alt$REF == ".")
d = which(alt$ALT == ".")
del = sort(unique(c(a,b,c,d)))
if (length(del) > 0) {
	alt = alt[-del, ]
}
cat("OK\n")



cat("Adjusting genotypes ... ")
del = unique(c(grep("^\\.[\\/\\|].+$", ref$SAMPLE), grep("^[01\\.][\\/\\|]\\..*$", ref$SAMPLE)))
if (length(del) > 0) {
	ref = ref[-del, ]
}

del = unique(c(grep("^\\.[\\/\\|].+$", alt$SAMPLE), grep("^[01\\.][\\/\\|]\\..*$", alt$SAMPLE), grep("^[01\\.]:$", alt$SAMPLE)))
if (length(del) > 0) {
	alt = alt[-del, ]
}

ref$GT = as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\1", ref$SAMPLE)) + as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\2", ref$SAMPLE))
alt$GT = as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\1", alt$SAMPLE)) + as.numeric(sub("^([01])[\\/\\|]([01]).*$", "\\2", alt$SAMPLE))
cat("OK\n")



cat("Generating keys ... ")
ref$KEY = paste(ref$CHROM, as.character(ref$POS), sep = "_")
alt$KEY = paste(alt$CHROM, as.character(alt$POS), sep = "_")

del = which(duplicated(ref$KEY))
if (length(del) > 0) {
	ref = ref[-del, ]
}

del = which(duplicated(alt$KEY))
if (length(del) > 0) {
	alt = alt[-del, ]
}
cat("OK\n")



key = intersect(ref$KEY, alt$KEY)
cat("Number of corresponding variants:", length(key), "\n")



cat("Subsetting data ... ")
sub.ref = ref[which(ref$KEY %in% key), ]
sub.alt = alt[which(alt$KEY %in% key), ]

sub.ref = sub.ref[order(sub.ref$KEY), ]
sub.alt = sub.alt[order(sub.alt$KEY), ]

if (!identical(sub.ref$POS, sub.alt$POS)) {
	stop("Wrong order!")
}
cat("OK\n")



cat("Setting hom.ref ")
hom.alt = alt[which(!alt$KEY %in% key), ]
hom.ref = NULL

chr.bed = split(bed, bed$CHROM)
chr.alt = split(hom.alt, hom.alt$CHROM)

for (chr in names(chr.alt)) {
	cat(".")
	ca = chr.alt[[chr]]
	cb = chr.bed[[chr]]
	
	x = rep(FALSE, cb$END[nrow(cb)])
	for (i in 1:nrow(cb)) {
		x[(cb$BEG[i]):(cb$END[i])] = TRUE
	}
	y = which(x[ca$POS])
	
	hom.ref = rbind(hom.ref, ca[y, ])
}

hom.ref$GT = 0
hom.alt = hom.alt[which(hom.alt$KEY %in% hom.ref$KEY), ]

hom.ref = hom.ref[order(hom.ref$KEY), ]
hom.alt = hom.alt[order(hom.alt$KEY), ]

if (!identical(hom.ref$POS, hom.alt$POS)) {
	stop("Wrong order!")
}
cat(" OK\n")



Ref = rbind(sub.ref, hom.ref)
Alt = rbind(sub.alt, hom.alt)


del = which(Ref$REF != Alt$REF)
cat("Removing REF inconsistencies:", length(del), "\n")
if (length(del) > 0) {
	Ref = Ref[-del, ]
	Alt = Alt[-del, ]
}

del = which(Ref$ALT != Alt$ALT)
cat("Removing ALT inconsistencies:", length(del), "\n")
if (length(del) > 0) {
	Ref = Ref[-del, ]
	Alt = Alt[-del, ]
}



cat("Compiling error profile ... \n")

p = data.frame(chr = Ref$CHROM, 
							 pos = Ref$POS, 
							 a0 = Ref$REF, 
							 a1 = Ref$ALT, 
							 maf = NA, 
							 ref.gt = Ref$GT,
							 alt.gt = Alt$GT,
							 error = (Ref$GT != Alt$GT),
							 type = NA,
							 label = NA,
							 order = NA,
							 stringsAsFactors = F)
rownames(p) = Ref$KEY

x = intersect(rownames(p), names(maf))

cat("Inserting MAF at n sites:", length(x), "\n")

if (!identical(rownames(p), x)) {
	p = p[x, ]
}

p$maf = maf[x]


del = which(p$maf == 0 & !p$error)
if (length(del) > 0) {
	p = p[-del, ]
}


p$type[ which(p$ref.gt == 0 & p$alt.gt == 0) ] = "00" # correct
p$type[ which(p$ref.gt == 0 & p$alt.gt == 1) ] = "01"
p$type[ which(p$ref.gt == 0 & p$alt.gt == 2) ] = "02"

p$type[ which(p$ref.gt == 1 & p$alt.gt == 0) ] = "10"
p$type[ which(p$ref.gt == 1 & p$alt.gt == 1) ] = "11" # correct
p$type[ which(p$ref.gt == 1 & p$alt.gt == 2) ] = "12"

p$type[ which(p$ref.gt == 2 & p$alt.gt == 0) ] = "20"
p$type[ which(p$ref.gt == 2 & p$alt.gt == 1) ] = "21"
p$type[ which(p$ref.gt == 2 & p$alt.gt == 2) ] = "22" # correct


p$label[ which(p$ref.gt == 0 & p$alt.gt == 0) ] = "Truth 0, called 0" # correct
p$label[ which(p$ref.gt == 0 & p$alt.gt == 1) ] = "Truth 0, called 1"
p$label[ which(p$ref.gt == 0 & p$alt.gt == 2) ] = "Truth 0, called 2"

p$label[ which(p$ref.gt == 1 & p$alt.gt == 0) ] = "Truth 1, called 0"
p$label[ which(p$ref.gt == 1 & p$alt.gt == 1) ] = "Truth 1, called 1" # correct
p$label[ which(p$ref.gt == 1 & p$alt.gt == 2) ] = "Truth 1, called 2"

p$label[ which(p$ref.gt == 2 & p$alt.gt == 0) ] = "Truth 2, called 0"
p$label[ which(p$ref.gt == 2 & p$alt.gt == 1) ] = "Truth 2, called 1"
p$label[ which(p$ref.gt == 2 & p$alt.gt == 2) ] = "Truth 2, called 2" # correct


p$order[ which(p$ref.gt == 0 & p$alt.gt == 0) ] = "MATCH: Truth 0, called 0" # correct
p$order[ which(p$ref.gt == 0 & p$alt.gt == 1) ] = "ERROR: Truth 0, called 1"
p$order[ which(p$ref.gt == 0 & p$alt.gt == 2) ] = "ERROR: Truth 0, called 2"

p$order[ which(p$ref.gt == 1 & p$alt.gt == 0) ] = "ERROR: Truth 1, called 0"
p$order[ which(p$ref.gt == 1 & p$alt.gt == 1) ] = "MATCH: Truth 1, called 1" # correct
p$order[ which(p$ref.gt == 1 & p$alt.gt == 2) ] = "ERROR: Truth 1, called 2"

p$order[ which(p$ref.gt == 2 & p$alt.gt == 0) ] = "ERROR: Truth 2, called 0"
p$order[ which(p$ref.gt == 2 & p$alt.gt == 1) ] = "ERROR: Truth 2, called 1"
p$order[ which(p$ref.gt == 2 & p$alt.gt == 2) ] = "MATCH: Truth 2, called 2" # correct


save(p, file = out.file)

cat("OK\n")




