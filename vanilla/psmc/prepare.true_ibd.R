
library(data.table)


args = commandArgs(T)

file = args[1]

ccf = fread(file, header = T, stringsAsFactors = F)


cols = c("MarkerID", "Fk", "SampleID0", "Chr0", "SampleID1", "Chr1", "Shared")


d = ccf[, cols, with=F]
d = unique(d)

write.table(d, file = sprintf("_%s", basename(file)), append = F, quote = F, row.names = F, col.names = T)


