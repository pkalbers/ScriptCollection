
args = commandArgs(T)

map = read.table(args[1], header = T, stringsAsFactors = F)

map = map[, c(2:4)]

names(map) = c("pposition", "rrate", "gposition")

write.table(map, file = sprintf("shapeit.3col.%s", args[1]), quote = F, col.names = T, row.names = F)


