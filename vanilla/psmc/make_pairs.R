
library(data.table)


args = commandArgs(T)

file = args[1]

ccf = fread(file, header = T, stringsAsFactors = F)

M = fread("../vanilla.marker.txt", header = T, stringsAsFactors = F)

pairs = data.table(mid = ccf$MarkerID, 
									 pos = M$Position[ ccf$MarkerID + 1 ],
									 id0 = (ccf$SampleID0 * 2) + ccf$Chr0,
									 id1 = (ccf$SampleID1 * 2) + ccf$Chr1,
									 shr = ccf$Shared)

pairs = unique(pairs)


write.table(pairs, file = sprintf("pairs_%s", file), append = F, sep = " ", quote = F, row.names = F, col.names = T)


