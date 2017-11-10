
library(data.table)


files = dir(pattern = ".+\\.marker\\.txt$")

m = list()
for (file in files) {
	print(file)
	m[[file]] = fread(file)
}

for (i in 1:length(m)) {
	m[[i]]$key = sprintf("%d %d %d", m[[i]]$MarkerID, m[[i]]$AlleleCount1, m[[i]]$GenotypeCount1)
}

key = m[[1]]$key

for (i in 2:length(m)) {
	key = intersect(key, m[[i]]$key)
}

for (i in 1:length(m)) {
	x = which(m[[i]]$key %in% key)
	m[[i]] = m[[i]][x, ]
}


span = 2:25
nmax = 10000


marker = m[[1]]


x = which(marker$AlleleCount1 %in% span & marker$GenotypeCount1 %in% span)

pos = marker$Position[x]

if (length(pos) > nmax) {
	pos = sample(pos, nmax)
}

cat(pos, file = sprintf("expos.%d.txt", nmax), sep = "\n", append = F)


POS = split(pos, cut(1:length(pos), 100))


for (i in 1:length(POS)) {
	cat(POS[[i]], file = sprintf("./packs_hmm/expos.%d.pack%03d.txt", nmax, i), sep = "\n", append = F)
}



