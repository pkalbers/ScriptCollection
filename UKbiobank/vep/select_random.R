


library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	m = fread(sprintf("ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	m = m[which(m$AlleleCount1 > 1 & m$AlleleCount1 < 1515), ]
	m = m[which(m$AlleleCount1 == m$GenotypeCount1), ]
	
	t = table(m$AlleleCount1)
	p = exp(t[as.character(m$AlleleCount1)] / max(t))
	
	m = m[sample(1:nrow(m), size = 10000, replace = F, prob = p),]
	
	cat(sample(m$Position), sep = "\n", file = sprintf("selected_chr%d.random.txt", i))
	
}



for (i in 1:22) {
	
	x = readLines(sprintf("selected_chr%d.random.txt", i))
	
	
	m = fread(sprintf("ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	m = m[which(m$AlleleCount1 > 1 & m$AlleleCount1 < 1515), ]
	m = m[which(m$AlleleCount1 == m$GenotypeCount1), ]
	
	y = sample(m$Position, size = 1000, replace = F)
	
	pos = unique(c(x, y))
	pos = sample(pos)
	cat(pos, sep = "\n", file = sprintf("selected_chr%d.random.txt", i), append = F)
}




x = d[which(d$Impact == "PATHOGENIC"),]
x = split(x, x$Chr)

for (chr in x) {
	cat(chr$Chr[1], "\n")
	i = chr$Chr[1]
	pos = as.numeric(readLines(sprintf("selected_chr%d.pathog.txt", i)))
	
	m = fread(sprintf("ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	z = match(chr$MarkerID, m$MarkerID)
	done = m$Position[z]
	
	if (any(pos %in% done)) {
		
		del = which(pos %in% done)
		pos = pos[-del]
		pos = unique(pos)
		cat(pos, sep = "\n", file = sprintf("selected_chr%d.pathog.txt", i), append = F)
	}
}

