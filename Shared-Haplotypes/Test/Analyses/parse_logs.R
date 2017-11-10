
tab <- c()

for (file in dir("./logs/",pattern="^log", full.names=TRUE)) {
	chr <- as.numeric(sub(".+\\.chr([0-9]+)\\..+", "\\1", file))
	
	rare <- NA
	nvar <- NA
	scan <- NA
	time <- NA
	
	lines <- readLines(file)
	
	for (line in lines) {
		hash <- substr(line, 1, 1)
		if (hash != "#" && hash != "C")
			next
		
		tag <- substr(line, 3, 6)
		
		if (tag == "Vari") {
			nvar <- as.numeric(sub("[^0-9]+ ([0-9]+)", "\\1", line))
			next
		}
		
		if (tag == "Rare") {
			rare <- as.numeric(sub("[^0-9]+ ([0-9]+)", "\\1", line))
			next
		}
		
		if (tag == "Scan") {
			scan <- as.numeric(sub("[^0-9]+ ([0-9]+)", "\\1", line))
			next
		}
		
		if (tag == "U  t") {
			time <- as.numeric(sub("[^0-9]+ ([0-9]+) .+", "\\1", line))
			next
		}
	}
	
	tab <- rbind(tab, c(chr=chr, nvar=nvar, rare=rare, scan=scan, time=time))
}

tab <- tab[order(tab[, "chr"]), ]
tab[, "time"] <- tab[, "time"] / 60


