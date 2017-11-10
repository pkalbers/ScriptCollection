#
# read phased haplotypes
#

library(data.table)

args = commandArgs(T)

data.file = args[1]
phap.file = args[2]


load(data.file)


# get phased haplotypes
cat("Reading phased haplotypes\n")
tmp = readLines(phap.file, n = -1)
tmp = strsplit(tmp, " ", T)
tmp = lapply(tmp, function(x) as.integer(x[6:length(x)]))
PH = matrix(0, nrow = length(tmp), ncol = length(tmp[[1]]))
for (i in 1:length(tmp)) {
	PH[i, ] = tmp[[i]]
}


