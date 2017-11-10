
library(data.table)

args = commandArgs(TRUE)
file = args[1]

#data = read.table(file, header=TRUE, stringsAsFactors=FALSE)
data = fread(file, header=TRUE, stringsAsFactors=FALSE)

save(data, file = sprintf("data.%s.RData", file))


