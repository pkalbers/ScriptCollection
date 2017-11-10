
library(data.table)

args = commandArgs(T)

d = fread(args[1], header = F, stringsAsFactors = F)

print(range(d$V6))
print(range(d$V7))

print(fivenum(d$V7 - d$V6 + 1))

key = sprintf("%s %s %d %d", d$V1, d$V3, d$V6, d$V7)
del = which(duplicated(key))

cat("N:", length(key), "dup:", length(del), "\n")

d = d[-del,]

print(fivenum(d$V7 - d$V6 + 1))


