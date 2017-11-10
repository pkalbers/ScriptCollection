

library(data.table)


nex = 5
mex = 1000


load("../mask/pos.chr20.filtered.RData")


x = which(p$dp >= quantile(p$dp, probs = 0.95))
if (length(x) == 0) stop("NONE!")
p = p[x,]


x = which(p$ns >= quantile(p$ns, probs = 0.95))
if (length(x) == 0) stop("NONE!")
p = p[x,]


x = which(p$AlleleCount1 >= 2 & p$AlleleCount1 <= max(p$AlleleCount1) / 2)
if (length(x) == 0) stop("NONE!")
p = p[x,]


#q = p[sample(1:nrow(p), size = mex, replace = F), ]
q = p[sample(1:nrow(p), size = mex, replace = F, prob = sqrt(p$AlleleCount1)), ]
#q = p[sample(1:nrow(p), size = mex, replace = F, prob = p$AlleleCount1), ]

plot(table(q$AlleleCount1))


#cat(sort(q$Position), file = "chr20_expos.1000.txt", sep = "\n")

pos = q$Position


tmp = split(pos, cut(1:length(pos), breaks = mex/nex))

na = lapply(tmp, function(pos) {
	if (length(pos) > 0) {
		cat(sort(pos), file = sprintf("expos/pack.%s.txt", paste(sample(c(0:9, LETTERS, letters), 16), collapse = "")), sep = "\n")
	}
})






