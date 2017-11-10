#
# shared haplotype breakpoint detection
#

library(data.table)

fk.max = 25

hap.file = "sim2000.hap.byrow.hap"
pos.file = "sim2000.pos.pos"



# get simulated haplotypes

tmp = readLines(hap.file, n = -1)
tmp = strsplit(tmp, "", T)
tmp = lapply(tmp, as.integer)
H = matrix(0, nrow = length(tmp[[1]]), ncol = length(tmp))
for (i in 1:length(tmp)) {
	H[, i] = tmp[[i]]
}

# make genotypes

G = matrix(0, nrow = nrow(H), ncol = ncol(H) / 2)
j = 1
for (i in 1:ncol(G)) {
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
}

# get positions

P = as.numeric(readLines(pos.file, n = -1))

# get frequencies

mac = rowSums(H)
maf = mac / ncol(H)
i = which(mac > (ncol(H) / 2))
if (length(i) > 0) {
	mac[i] = ncol(H) - mac[i]
	maf[i] = 1 - maf[i]
}

# make fk index

fki = which(mac > 1 & mac <= fk.max)
fki = lapply(fki, function(i) {
	s = which(G[i, ] == 1)
	n = length(s)
	if (n == 1) return(NULL)
	data.table(index = i, 
						 position = P[i], 
						 n.sharer = n, 
						 sharer = paste(as.character(s), collapse = "|"),
						 is.double.het = (mac[i] != n))
})
fki = rbindlist(fki)

# function to parse identified sharers

parse.sharers = function(x) {
	as.integer(strsplit(x, "|", T)[[1]])
}

# function to address an individual's haplotypes in matrix

get.indv.haps = function(i) {
	c(i * 2 - 1, i * 2)
}

# save data

save(H, G, P, 
		 mac, maf, fki, 
		 parse.sharers,
		 get.indv.haps,
		 file = "sim2000.RData")




#
# write sharing pair table for IBD detection
#

ibd.file = file("sim2000.pairs.txt", open = "w")

count = 0
for (i in 1:nrow(fki)) {
	idx = fki$index[i]
	pos = fki$position[i]
	set = parse.sharers(fki$sharer[i])
	set = combn(set, 2, simplify = F)
	for (pair in set) {
		cat(sprintf("%d %.8f %d %d\n", idx, pos, pair[1], pair[2]), file = ibd.file)
		count = count + 1
	}
}

cat("Number of pairs: ", count, "\n")

close(ibd.file)



#
# write GEN + sample file for phasing in shapeit
#













