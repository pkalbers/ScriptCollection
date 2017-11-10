
library(data.table)

args = commandArgs(T)

file = args[1]
pack = as.numeric(args[2])

pairs = fread(file, header = T)

size = as.integer(round(nrow(pairs) / pack))

dir.create("packs")

a = 1
b = size

for (i in 1:pack) {
	tmp = NULL
	print(c(a, b))

	if (i == pack) {
		tmp = pairs[a:nrow(pairs), ]
	} else {
		tmp = pairs[a:b, ]
	}

	if (nrow(tmp) == 0) break

	write.table(tmp, file = sprintf("./packs/pack.%04d.%s", i, file), append = F, quote = F, row.names = F, col.names = T)

	a = b + 1
	b = a + size
}


