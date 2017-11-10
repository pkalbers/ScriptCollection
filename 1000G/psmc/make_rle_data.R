

library(data.table)


size = 100000
skip = 0
i = 0
repeat {
	i = i + 1
	
	cat(i, "...")
	tmp = fread("1000G_chr20.impute.hap", sep = " ", header = F, stringsAsFactors = F, colClasses = "integer", nrows = size, skip = skip)
	cat(" OK\n")
	
	save(tmp, file = sprintf("tmp%05d.RData", i))
	
	if (nrow(tmp) < size) {
		cat(" DONE")
		break
	}
	
	skip = skip + size
}




files = sort(dir(pattern = "^tmp([0-9]+)\\.RData$"))

for (file in files) {
	cat(file, "\n")
	
	load(file)
	
	tmp = as.list(tmp)
	tmp = lapply(tmp, rle)
	
	save(tmp, file = sprintf("rle_%s", file))
}




files = sort(dir(pattern = "^rle_tmp([0-9]+)\\.RData$"))

d = list()

for (file in files) {
	cat(file, "\n")
	
	load(file)
	
	d = c(d, list(tmp))
}



num = unique(sapply(d, length))
tag = names(d[[1]])


q = list()

for (i in 1:num) {
	cat(".")
	
	h = c()
	
	for (j in 1:length(d)) {
		r = d[[j]][[tag[i]]]
		r = inverse.rle(r)
		h = c(h, r)
	}
	
	q[[tag[i]]] = rle(h)
}
cat("\n")



# mxl = max(sapply(q, function(x) length(x$lengths)))
# 
# 
# L = matrix(NA, nrow = mxl, ncol = num)
# V = matrix(NA, nrow = mxl, ncol = num)
# 
# 
# for (i in 1:num) {
# 	if (i %% 500 == 0) cat(i, "\n")
# 	
# 	p = q[[tag[i]]]
# 	x = length(p$lengths)
# 	
# 	L[1:x, i] = p$lengths
# 	V[1:x, i] = p$values
# }



H = q
rm(q)
gc()



a = fread("1000G_chr20.impute.legend", header = T, stringsAsFactors = F)
b = fread("../data/1000G.chr20.marker.txt", header = T, stringsAsFactors = F)



z = which(a$pos %in% b$Position)
z = intersect(z, which(nchar(a$allele0) == 1))
z = intersect(z, which(nchar(a$allele1) == 1))

table(a$allele0[z])
table(a$allele1[z])

all(a$pos[z] %in% b$Position)
identical(a$pos[z], b$Position)

pos = b$Position

z = sort(z)

H = lapply(H, function(h, k) {
	rle(inverse.rle(h)[k])
}, z)


save(H, num, tag, pos, file = "snp_1000G_chr20_rle_data.RData")















