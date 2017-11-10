
args = commandArgs(TRUE)

#
# load mask
#

mask_file = args[1]

mask = read.table(mask_file, header = FALSE)

mask = split(mask, mask$V1)

mask = mask$chr20 ##########
mask$V1 <- NULL
mask$V4 <- NULL
rownames(mask) <- NULL
names(mask) <- c("from", "to")


make.is.masked <- function(m) {
	mask <- m
	f <- function(pos) {
		return(length(which(mask$from <= pos & mask$to >= pos)) != 0)
	}
	return( f )
}

is.masked <- make.is.masked(mask)



#
# filter markers using mask
#

marker_file = args[2]

marker = read.table(marker_file, header = TRUE, stringsAsFactors = FALSE)

masked = rep(NA, nrow(marker))
for (i in 1:nrow(marker)) {
	masked[i] = is.masked(marker$position[i])
}

# write sites to exclude
excl = marker$position[!masked]

writeLines(as.character(excl), con = "../excluded_sites.txt")

#save(is.masked, file = "mask_function.RData")
