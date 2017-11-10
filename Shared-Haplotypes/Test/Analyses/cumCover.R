
load("allChr.shap.RData")

library("parallel")

seq.from <- min(shap$from)
seq.to <- max(shap$to)

splt <- split(shap, shap$mac)

from <- mclapply(splt, function(s) split(s, s$from), mc.cores=2)

maxdist < -mclapply(from, function(f) lapply(f, function(t) max(t$to)), mc.cores=2)

save(maxdist, file="plot.maxdist.raw.RData")

M <- c()
for (mac.tag in names(maxdist)) {
	mac <- as.numeric(mac.tag)
	
	for (fr.tag in names(maxdist[[mac.tag]])) {
		fr <- as.numeric(fr.tag)
		to <- maxdist[[mac.tag]][[fr.tag]]
		
		M <- rbind(M, c(mac=mac, from=fr, to=to))
	}
}

save(M, file="plot.maxdist.res.RData")
