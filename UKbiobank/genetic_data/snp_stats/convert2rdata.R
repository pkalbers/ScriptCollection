#
# convert stats to RData
#


stat.files = dir(pattern = "snp_stats_chr[0-9]+.txt.gz")

for (stat.file in stat.files) {
	
	cat(stat.file, "\n")
	
	stat = read.table(stat.file, header = T, stringsAsFactors = F)
	
	save(stat, file = sub("^(.+)txt\\.gz$", "\\1RData", stat.file))
	
}










