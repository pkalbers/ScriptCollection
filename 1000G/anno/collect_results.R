
library(data.table)


files = dir(pattern = "^result\\.([A-Z]{3})\\.chr([0-9]+)\\..+\\.age\\.sites\\..+\\.txt$", full.names = T, recursive = T)

d = NULL

for (file in files) {
	str = basename(file)

	cat(str, " ")

	tmp = fread(file, header = T)

	if (nrow(tmp) == 0) {
		cat("\n")
		next
	}

	chr  = as.numeric(sub("^result\\.([A-Z]{3})\\.chr([0-9]+)\\.([0-9a-zA-Z_]+)\\.([0-9]+)\\.txt\\.age\\.sites\\..+\\.txt$", "\\2", str))
	type = sub("^result\\.([A-Z]{3})\\.chr([0-9]+)\\.([0-9a-zA-Z_]+)\\.([0-9]+)\\.txt\\.age\\.sites\\..+\\.txt$", "\\1", str)
	mode = sub("^result\\.([A-Z]{3})\\.chr([0-9]+)\\.([0-9a-zA-Z_]+)\\.([0-9]+)\\.txt\\.age\\.sites\\..+\\.txt$", "\\3", str)

	cat(type, "")
	cat(mode, "")

	tmp$chr = chr
	tmp$type = type
	tmp$mode = mode

	d = rbind(d, tmp)

	cat("\n")
}

save(d, file = "results.collected.RData")

