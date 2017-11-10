#
# apply error profile
#

args = commandArgs(T)

data.file = args[1]
errp.file = args[2] ### error matrix saved in 'm'
errp.name = args[3]

load(data.file)
load(errp.file)

prefix = sub("^(.+)\\.RData$", "\\1", data.file)


# convert genotyping error to haplotypes

err = matrix(0, 4, 4, dimnames = list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))

err["00", "00"] = m["0", "0"]
err["00", "01"] = m["0", "1"] / 2
err["00", "10"] = m["0", "1"] / 2
err["00", "11"] = m["0", "2"]

err["01", "00"] = m["1", "0"]
err["01", "01"] = m["1", "1"]
err["01", "10"] = 0           # otherwise read as flip error
err["01", "11"] = m["1", "2"]

err["10", "00"] = m["1", "0"]
err["10", "01"] = 0           # otherwise read as flip error
err["10", "10"] = m["1", "1"]
err["10", "11"] = m["1", "2"]

err["11", "00"] = m["2", "0"]
err["11", "01"] = m["2", "1"] / 2
err["11", "10"] = m["2", "1"] / 2
err["11", "11"] = m["2", "2"]


### apply error to haplotypes
cat("Apply error to haplotypes:\n")
j = 1
for (i in 1:ncol(G)) {
	if (i %% 100 == 0) cat(".")
	
	hm = H[, j]
	hp = H[, j+1]

	h = paste(as.character(hm), as.character(hp), sep="")
	
	h00 = which(h == "00")
	h01 = which(h == "01")
	h10 = which(h == "10")
	h11 = which(h == "11")
	
	h[h00] = sample(c("00", "01", "10" ,"11"), size = length(h00), replace = T, prob = err["00", ])
	h[h01] = sample(c("00", "01", "10" ,"11"), size = length(h01), replace = T, prob = err["01", ])
	h[h10] = sample(c("00", "01", "10" ,"11"), size = length(h10), replace = T, prob = err["10", ])
	h[h11] = sample(c("00", "01", "10" ,"11"), size = length(h11), replace = T, prob = err["11", ])
	
	hm = as.numeric(substr(h, 1, 1))
	hp = as.numeric(substr(h, 2, 2))
	
	H[, j]   = hm
	H[, j+1] = hp
	
	j = j + 2
}
cat("\n")



# make genotypes
# keep track of introduced error
cat("Make genotypes from errorised haplotypes:\n")
E = m
E[, ] = 0
j = 1
for (i in 1:ncol(G)) {
	if (i %% 100 == 0) cat(".")
	
	g0 = which(G[, i] == 0)
	g1 = which(G[, i] == 1)
	g2 = which(G[, i] == 2)
	
	G[, i] = H[, j] + H[, j+1]
	j = j + 2
	
	E["0", "0"] = E["0", "0"] + sum(G[g0, i] == 0)
	E["0", "1"] = E["0", "1"] + sum(G[g0, i] == 1)
	E["0", "2"] = E["0", "2"] + sum(G[g0, i] == 2)
	
	E["1", "0"] = E["1", "0"] + sum(G[g1, i] == 0)
	E["1", "1"] = E["1", "1"] + sum(G[g1, i] == 1)
	E["1", "2"] = E["1", "2"] + sum(G[g1, i] == 2)
	
	E["2", "0"] = E["2", "0"] + sum(G[g2, i] == 0)
	E["2", "1"] = E["2", "1"] + sum(G[g2, i] == 1)
	E["2", "2"] = E["2", "2"] + sum(G[g2, i] == 2)
}
E = E / rowSums(E)
cat("\n")


# get frequencies
mac = rowSums(H)
maf = mac / ncol(H)
frq = maf
i = which(mac > (ncol(H) / 2))
if (length(i) > 0) {
	mac[i] = ncol(H) - mac[i]
	maf[i] = 1 - maf[i]
}


# make fk index
fki = which(mac > 1 & mac <= fk.max)
fki = lapply(fki, function(i) {
	g = which(G[i, ] == 1)
	h = if (frq[i] <= 0.5) which(H[i, ] == 1) else which(H[i, ] == 0)
	h = intersect(h, get.indv.haps(g))
	n = length(g)
	if (n < 2) return(NULL)
	data.table(index = i, 
						 pos = P[i], 
						 frq = frq[i],
						 maf = maf[i],
						 mac = mac[i],
						 n   = n, 
						 g.sharer = paste(as.character(g), collapse = "|"),
						 h.sharer = paste(as.character(h), collapse = "|"))
})
fki = rbindlist(fki)


# save data
save(H, G, P, E, 
		 mac, maf, fki, 
		 parse.sharers,
		 get.indv.haps,
		 file = sprintf("%s.%s.RData", prefix, errp.name))




