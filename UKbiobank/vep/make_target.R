
args = commandArgs(T)

n = as.numeric(args[1])


file = sample(dir(pattern = "^selected_chr.+"))[1]

chr = as.numeric(sub("^selected_chr([0-9]+)\\..+", "\\1", file))
imp = sub(".+impact\\.([A-Z]+)\\.txt$", "\\1", file)


pos = readLines(file)
pos = sample(pos, size = n)


str = paste(sample(c(0:9, letters, LETTERS), 16, replace=TRUE), collapse="") 

out = sprintf("target.chr%d.%s.%s.txt", chr, imp, str)

cat(pos, sep = "\n", file = out)

cat(chr, imp, out, str)


