#
# write sharing pair table for Shared Haplotype (approx IBD) detection
#

library(data.table)

#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20.GenErr_1000G" #args[1]
chunks = 250 # args[2]

load(sprintf("data.%s.RData", prefix))


dir.create("pairs", showWarnings = F)
setwd("pairs")


cat("Distributing sharing pair tables ... \n")

splt = split(FKI, cut(1:nrow(FKI), chunks))


for (i in 1:length(splt)) {
	cat(sprintf(" %d of %d", i, length(splt)))
	
	pair.name = sprintf("pairs.%s.%04d.RData", prefix, i)
	
	pair = splt[[i]]
	
	save(pair, file = pair.name)
	cat("\n")
}

cat("Done\n")

