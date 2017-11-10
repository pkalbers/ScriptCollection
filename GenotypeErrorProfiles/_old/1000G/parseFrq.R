#
# parse frequencies
#

library(data.table)


frq.files = dir(pattern = ".+frq$")

frq = list()

for (frq.file in frq.files) {
	
	cat("Reading:", frq.file, "... ")
	tmp = readLines(frq.file)
	cat("OK\n")
	
	tmp = tmp[-1]
	
	
	cat("Splitting strings ... ")
	tmp = strsplit(tmp, "\t")
	cat("OK\n")
	
	
	cat("Generating keys ... ")
	key = lapply(tmp, function(x) {
		paste(x[1], x[2], sep = "_")
	})
	key = unlist(key)
	cat("OK\n")
	
	names(tmp) = key
	
	
	cat("Parsing frequencies ... ")
	str = sapply(tmp, function(x) {
		x = x[-(1:4)]
		if (length(x) > 2) return(NULL)
		paste(x, collapse = "__")
	})
	
	
	tab = data.frame(key = names(str), 
									 chr = as.numeric(sub(".+\\.chr([0-9]+)\\..+", "\\1", frq.file)), 
									 pos = as.numeric(sub("^[0-9]+_([0-9]+)$", "\\1", names(str))), 
									 ref = sub("([ACGT]):([0-9\\.]+)__([ACGT]):([0-9\\.]+)", "\\1", str),
									 alt = sub("([ACGT]):([0-9\\.]+)__([ACGT]):([0-9\\.]+)", "\\3", str),
									 ref.frq = as.numeric(sub("([ACGT]):([0-9\\.]+)__([ACGT]):([0-9\\.]+)", "\\2", str)),
									 alt.frq = as.numeric(sub("([ACGT]):([0-9\\.]+)__([ACGT]):([0-9\\.]+)", "\\4", str)),
									 stringsAsFactors = F)
	
	a = which(is.na(tab$chr))
	b = which(is.na(tab$pos))
	c = which(is.na(tab$ref.frq))
	d = which(is.na(tab$ref.frq))
	e = which(! tab$ref %in% c("A", "C", "G", "T"))
	f = which(! tab$alt %in% c("A", "C", "G", "T"))
	
	del = unique(c(a, b, c, d, e, f))
	if (length(del) > 0) {
		tab = tab[-del, ]
	}
	
	del = which(duplicated(tab$key))
	if (length(del) > 0) {
		tab = tab[-del, ]
	}
	
	tab$alt.is.minor = T
	tab$maf = tab$alt.frq
	x = which(tab$ref.frq < tab$alt.frq)
	if (length(x) > 0) {
		tab$alt.is.minor[x] = F
		tab$maf[x] = tab$ref.frq[x]
	}
	
	frq[[frq.file]] = as.data.table(tab)
	
	cat("OK\n\n")
}


frq = rbindlist(frq)
frq = as.data.frame(frq)
rownames(frq) = frq$key

cat("Saving data ... ")
save(frq, file = "frq.RData")
cat("DONE\n")











