

args <- commandArgs(TRUE)
vcf <- args[1]
nlines <- args[2]

hapsplit <- function(tok) 
{ 
	a <- substr(tok, 1, 1)
	b <- substr(tok, 3, 3)
	c(a, b)
}

con <- gzcon(file(vcf, "r"))

raw <- list()
n <- 1
count <- 0
while(TRUE)
{
	lines <- readLines(con, n)
	
	if (length(lines) == 0) break
	
	if (n == 1)
	{
		if (substr(lines, 1, 1) == "#") next else n <- nlines
	}
	
	lines <- lapply(lines, function(line) unlist(strsplit(line, "\t", TRUE))) #, mc.cores = length(lines))
	
	for (i in 1:length(lines))
	{
		if (lines[[i]][7] != "PASS") next
		
		x <- table(unlist(lapply(lines[[i]][10:length(lines[[i]])], hapsplit))) #, mc.cores = length(lines)))
		
		if ("." %in% names(x)) next
		
		c <- as.character(length(x))
		
		if (c == "1") next
		
		if (c %in% names(raw))
			raw[[c]] <- c(raw[[c]], list(x))
		else
		{
			raw[[c]] <- list()
			raw[[c]] <- c(raw[[c]], list(x))
		}
	}

	count <- count + length(lines)
	if (count %% 1000 == 0) cat(count, "\n")
}

res <- list()
tab <- list()

for (c in names(raw))
{
	res[[c]] <- unlist(lapply(raw[[c]], function(x) { x[which.min(x)] }), use.names = FALSE)
	tab[[c]] <- table(res[[c]])
}

save(raw, res, tab, file = sprintf("_%s_.RData", vcf))


