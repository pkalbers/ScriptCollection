


parse.genotypes <- function(data) {
	cols = grep("^genotype.+", names(data), value = T)
	
	g = list()
	
	for (col in cols) {
		tags = strsplit(col, "_", T)[[1]]
		
		ind = tags[3]
		
		chr = data[[col]]
		
		x = which(chr == 10)
		if (length(x) > 0) chr[x] = 1
		
		x = which(chr == 11)
		if (length(x) > 0) chr[x] = 2
		
		g[[ind]] = chr
		
	}
	g
}


parse.target <- function(data) {
	col = grep("^target.+", names(data), value = T)
	
	tags = strsplit(col, "_", T)[[1]]
	
	ind = tags[3]
	
	chr = data[[col]]
	
	x = which(chr == 10)
	if (length(x) > 0) chr[x] = 1
	
	x = which(chr == 11)
	if (length(x) > 0) chr[x] = 2
	
	chr
}


parse.sample <- function(data) {
	cols = grep("^target|genotype.+", names(data), value = T)
	
	s = NULL
	
	for (col in cols) {
		tags = strsplit(col, "_", T)[[1]]
		
		key = tags[2]
		ind = tags[3]
		pop = tags[4]
		grp = tags[5]
		
		s = rbind(s, data.frame(key, ind, pop, grp))
	}
	rownames(s) = s$ind
	s
}


parse.sharedhap <- function(data) {
	cols = grep("^shared.+", names(data), value = T)
	
	h = list()
	
	for (col in cols) {
		tags = strsplit(col, "_", T)[[1]]
		
		ind = tags[3]
		
		sh = rep(FALSE, length(data[[col]]))
		
		x = which(data[[col]] != ".")
		if (length(x) > 0) sh[x] = TRUE
		
		h[[ind]] = sh
		
	}
	h
}


parse.meta <- function(data, shap) {
	meta = data[, 2:8]
	meta = cbind(meta, maf = apply(data.frame(data$frq0, data$frq1), 1, min), mac = apply(data.frame(data$count0, data$count1), 1, min), f = 0)
	
	shap = as.data.frame(shap)
	shap = apply(shap, 1, function(x) length(which(x)))
	shap = shap + 1
	
	meta$f = shap
	
	meta
}


args = commandArgs(TRUE)
file = args[1]

load(file)

indv = parse.sample(data)
shap = parse.sharedhap(data)
meta = parse.meta(data, shap)
target = parse.target(data)
others = parse.genotypes(data)

save(indv, shap, meta, target, others, file = sprintf("parsed.%s", file))

