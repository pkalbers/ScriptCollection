
library(data.table)


files = dir(pattern = "^count\\.chr[0-9]+\\.psm\\.[0-9]+\\.txt$")

x = dim(fread(files[1], header = T, stringsAsFactors = F))

d = list()

for (file in files) {
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	if (!identical(dim(tmp), x)) {
		print(file)
		next
	}
	
	chr = as.numeric(sub("^count\\.chr([0-9]+)\\.psm\\.([0-9]+)\\.txt$", "\\1", file))
	fk  = as.numeric(sub("^count\\.chr([0-9]+)\\.psm\\.([0-9]+)\\.txt$", "\\2", file))
	
	cat(sprintf(" %d  %d \n", chr, fk))
	
	if (!as.character(chr) %in% names(d)) {
		d[[as.character(chr)]] = list()
	}
	
	d[[as.character(chr)]][[as.character(fk)]] = tmp
}


save(d, file = "data.psm.RData")



#####



## merge by chromosome

load("data.psm.RData")

p = d[["1"]]

for (chr in as.character(2:22)) {
	cat(chr, " ")
	
	for (fk in names(p)) {
		cat(".")
		
		if (identical(p[[fk]]$SampleID0, d[[chr]][[fk]]$SampleID0) &&
				identical(p[[fk]]$SampleID1, d[[chr]][[fk]]$SampleID1)) {
			
			p[[fk]]$N = p[[fk]]$N + d[[chr]][[fk]]$N
			
			next
		}
		stop("!!!")
	}
	
	cat("\n")
}

d = p
save(d, file = "data.psm.byChr.RData")




#######


## merge by fk

load("data.psm.RData")



