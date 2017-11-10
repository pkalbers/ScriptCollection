

library(data.table)


nex = 5
mex = 1000


M = fread("../vanilla.marker.txt", header = T, stringsAsFactors = F)


x = which(M$AlleleCount1 >= 2 & M$AlleleCount1 <= max(M$AlleleCount1) / 2)
x = sample(x, size = mex, replace = F, prob = log(M$AlleleCount1[x]))


pos = sample(M$Position[x])

tmp = split(pos, cut(1:length(pos), breaks = mex/nex))

na = lapply(tmp, function(pos) {
	if (length(pos) > 0) {
		cat(sort(pos), file = sprintf("expos/pack.%s.txt", paste(sample(c(0:9, LETTERS, letters), 16), collapse = "")), sep = "\n")
	}
})






