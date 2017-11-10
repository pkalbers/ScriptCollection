


files <- dir(pattern = "runtime")

data <- data.frame(file = files, type = NA, f = NA, cpu = NA, user = NA, shs = NA, vm = NA, rss = NA, vm.gb = NA, rss.gb = NA)


for (i in 1:nrow(data))
{
	filename <- as.character(data[i, 1])
	
	l <- strsplit(filename, split = ".", fixed = TRUE)[[1]]
	
	data[i, 2] <- l[1]
	data[i, 3] <- as.numeric(l[2])
	
	line <- readLines(filename)
	
	data[i, 4] <- as.numeric(line[3])
	data[i, 5] <- as.numeric(line[2])
	
	data[i, 7] <- as.numeric(sub("^VM ([0-9]+)\\.[0-9]+$", "\\1", line[4]))
	data[i, 8] <- as.numeric(sub("^RSS ([0-9]+)\\.[0-9]+$", "\\1", line[5]))
	
	data[i, 9]  <- data[i, 7] / 1000000
	data[i, 10] <- data[i, 8] / 1000000
	
	lines <- readLines(sprintf("%s.%d.log", data[i, 2], data[i, 3]))
	
	for (line in lines)
	{
		if (grepl("^Identified", line))
		{
			data[i, 6] <- as.numeric(sub("^Identified rare haplotypes: ([0-9]+) \\(in [0-9]+ markers\\)$", "\\1", line))
		}
	}
}


splt <- split(data, data$type)

names(splt$cum) <- sprintf("cum.%s", names(splt$cum))
names(splt$exact) <- sprintf("exact.%s", names(splt$exact))

splt$cum$type <- NULL
splt$exact$type <- NULL


data <- cbind(splt$cum, splt$exact)


write.table(data, file = "plotdata.txt", quote = FALSE, row.names = FALSE)






######





files <- dir(pattern = "baseline.*memory")

data <- data.frame(file = files, type = NA, f = NA, vm = NA, rss = NA, vm.gb = NA, rss.gb = NA)

for (i in 1:nrow(data))
{
	filename <- as.character(data[i, 1])
	
	l <- strsplit(filename, split = ".", fixed = TRUE)[[1]]
	
	data[i, 2] <- l[2]
	data[i, 3] <- as.numeric(l[3])
	
	line <- readLines(filename)
	
	data[i, 4] <- as.numeric(sub("^VM ([0-9]+)\\.[0-9]+$", "\\1", line[1]))
	data[i, 5] <- as.numeric(sub("^RSS ([0-9]+)\\.[0-9]+$", "\\1", line[2]))
	
	data[i, 6] <- data[i, 4] / 1000000
	data[i, 7] <- data[i, 5] / 1000000
}





