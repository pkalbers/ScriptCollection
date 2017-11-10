#
# include Beagle IBD results
#

library(data.table)
library(Hmisc)


.prefix = "OutOfAfricaHapMap20.GenErr_1000G"
.gmap.file = "../genetic_map_GRCh37_chr20.txt"

load(sprintf("data.%s.RData", .prefix))


if ("GPOS" %in% ls()) {
	rm(GPOS)
}


.gmap = read.table(.gmap.file, header = T, stringsAsFactors = F)
names(.gmap) = c("chr", "pos", "rate", "dist")


RATE = approxExtrap(.gmap$pos, .gmap$rate, POS)$y


save(list = ls(), file = sprintf("data.%s.RData", .prefix))




