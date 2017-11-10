
library(data.table)


file = "vanilla.marker.txt"

span = 2:20
nmax = 10000


marker = fread(file, header = T)


x = which(marker$AlleleCount1 %in% span & marker$GenotypeCount1 %in% span)

pos = marker$Position[x]

if (length(pos) > nmax) {
	pos = sample(pos, nmax)
}

cat(pos, file = sprintf("expos.%d.txt", nmax), sep = "\n", append = F)


