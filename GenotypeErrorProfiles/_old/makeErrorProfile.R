#
# make GT error profile
#

args = commandArgs(T)

truth.file = args[1]
typed.file = args[2]



cat("Loading data ... ")
load(truth.file)
load(typed.file)
cat("OK\n")


bed = split(truth.bed, truth.bed$CHROM)



cat("Sites contained in both truth & typed data ... ")

key = intersect(typed$key, truth$key)
sub.typed = typed[key, ]
sub.truth = truth[key, ]
									
inv = which(sub.truth$ref == sub.typed$alt & sub.truth$alt == sub.typed$ref)
if (length(inv) > 0) {
	cat(sprintf("(inverted alleles: %d) ", length(inv)))
	sub.typed$gt = abs(sub.typed$gt - 2)
}

a = which(sub.truth$ref != sub.typed$ref)
b = which(sub.truth$alt != sub.typed$alt)
del = sort(unique(c(a, b)))
if (length(del) > 0) {
	cat(sprintf("(non-matching alleles: %d) ", length(del)))
	sub.typed = sub.typed[-del, ]
	sub.truth = sub.truth[-del, ]
}

error.match = data.frame(frq.ref = sub.truth$frq.ref, 
												 frq.alt = sub.truth$frq.alt, 
												 truth.gt = sub.truth$gt,
												 typed.gt = sub.typed$gt,
												 truth.is.confident = sub.truth$is.confident,
												 truth.is.assumed = F,
												 typed.is.assumed = F)

cat("OK\n")



cat("Sites contained only in typed data ")

error.typed.diff = NULL

key = setdiff(typed$key, truth$key)
if (length(key) > 0) {
	
	sub.typed = typed[key, ]
	
	error.typed.diff = data.frame(frq.ref = sub.typed$frq.ref, 
																frq.alt = sub.typed$frq.alt, 
																truth.gt = 0,
																typed.gt = sub.typed$gt,
																truth.is.confident = NA,
																truth.is.assumed = T,
																typed.is.assumed = F)
	
	
	
	tmp = strsplit(key, "_", T)
	tmp = lapply(tmp, as.integer)
	tmp = data.frame(KEY   = key, 
									 CHROM = sapply(tmp, function(x) x[1]), 
									 POS   = sapply(tmp, function(x) x[2]))
	tmp = split(tmp, tmp$CHROM)
	
	tags = intersect(names(bed), names(tmp))
	
	truth.mask = rep(FALSE, nrow(error.typed.diff))
	names(truth.mask) = key
	
	for (tag in tags) {
		cat(".")
		cv = tmp[[tag]]
		cb = bed[[tag]]
		
		end = max(max(cb$END), max(cv$POS))
		
		x = rep(FALSE, end)
		for (i in 1:nrow(cb)) {
			x[(cb$BEG[i]):(cb$END[i])] = TRUE
		}
		y = which(x[cv$POS])
		
		if (length(y) > 0) {
			truth.mask[ cv$KEY[y] ] = TRUE
		}
	}
	
	error.typed.diff$truth.is.confident = truth.mask
	
}

cat("OK\n")



cat("Sites contained only in truth data ... ")

error.truth.diff = NULL

key = setdiff(truth$key, typed$key)
if (length(key) > 0) {
	
	sub.truth = truth[key, ]
	
	error.truth.diff = data.frame(frq.ref = sub.truth$frq.ref, 
																frq.alt = sub.truth$frq.alt, 
																truth.gt = sub.truth$gt,
																typed.gt = 0,
																truth.is.confident = sub.truth$is.confident,
																truth.is.assumed = F,
																typed.is.assumed = T)
}

cat("OK\n")




cat("Saving ... ")

error = rbind(error.match, 
							error.typed.diff,
							error.truth.diff)

error$is.error = (error$truth.gt != error$typed.gt)

save(error, file = "error.RData")

cat("OK\n")

