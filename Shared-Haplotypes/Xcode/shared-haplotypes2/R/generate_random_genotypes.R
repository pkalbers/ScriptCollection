#
# generate random genotypes
#

args <- commandArgs(TRUE)
Nsamples <- as.numeric(args[1])
Nvariants <- as.numeric(args[2])
file <- args[3]

types <- c(" 0 0 1", " 0 1 0", " 1 0 0")
alleles <- c("A", "C", "G", "T")

cat("", file = file)

for (i in 1:Nvariants) {
	data <- paste(sample(types, Nsamples, replace = TRUE), collapse = "")
	info <- paste("---", sprintf("snp_%d", i), i, sample(alleles, 1), sample(alleles, 1), collapse = " ")
	
	cat(info, data, "\n", file = file, append = TRUE, sep = "")
	
	if (i %% 10000 == 0) cat(".")
}
