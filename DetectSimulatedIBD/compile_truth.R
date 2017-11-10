#
# compile truth
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)


#args = commandArgs(T)

prefix = "history"  #args[1]             # file prefix

load(sprintf("%s.RData", prefix))


index = 1:length(POS)
index.names = sprintf("%.5f", POS)
names(index) = index.names


wall.lhs = POS[10]
wall.rhs = POS[length(POS)-10]

tru.files = dir(pattern = "^ibd_mrca\\..+\\.txt$", recursive = T)

truth = list()

for (tru.file in tru.files) {
	tru = read.table(tru.file, header = F, col.names = c("index", "position", "fk", "h0", "h1", "lhs", "rhs"), skip = 1)
	tru = as.data.table(tru)
	
	tag = sub("^.+\\.([0-9]+)\\..+$", "\\1", tru.file)
	
	cat(tag, "\n")
	
	
	tru$position = POS[ tru$index ]
	
	idx.lhs = index[ (sprintf("%.5f", tru$lhs)) ]
	idx.rhs = index[ (sprintf("%.5f", tru$rhs)) ]
	
	nas.lhs = which(is.na(idx.lhs))
	nas.rhs = which(is.na(idx.rhs))
	
	if (length(nas.lhs) > 0) idx.lhs[nas.lhs] = sapply(tru$lhs[nas.lhs], function(x) which.min(abs(x - POS)) )
	if (length(nas.rhs) > 0) idx.rhs[nas.rhs] = sapply(tru$rhs[nas.rhs], function(x) which.min(abs(x - POS)) )
	
	tru$lhs = POS[ idx.lhs ]
	tru$rhs = POS[ idx.rhs ]
	
	# reaches end of chromosome
	
	tru$wall.lhs = F; wall = which(round(tru$lhs) < wall.lhs); if (length(wall) > 0) tru$wall.lhs[wall] = T
	tru$wall.rhs = F; wall = which(round(tru$rhs) > wall.rhs); if (length(wall) > 0) tru$wall.rhs[wall] = T
	
	
	# distance to focal site
	
	tru$dist.lhs = tru$lhs - tru$position
	tru$dist.rhs = tru$rhs - tru$position
	
	
	# breakpoint positions
	
	prebreak.lhs = sprintf("%.5f", tru$lhs[!tru$wall.lhs])
	prebreak.rhs = sprintf("%.5f", tru$rhs[!tru$wall.rhs])
	
	if (!all(prebreak.lhs %in% index.names) || !all(prebreak.rhs %in% index.names)) {
		stop("Cannot match all breakpoint positions")
	}
	
	tru$brk.lhs = NA;  tru$brk.lhs[!tru$wall.lhs] = POS[ index[prebreak.lhs] - 1 ]
	tru$brk.rhs = NA;  tru$brk.rhs[!tru$wall.rhs] = POS[ index[prebreak.rhs] + 1 ]
	
	
	# breakpoint distances
	
	tru$brk.dist.lhs = NA;  tru$brk.dist.lhs[!tru$wall.lhs] = tru$brk.lhs[!tru$wall.lhs] - tru$position[!tru$wall.lhs]
	tru$brk.dist.rhs = NA;  tru$brk.dist.rhs[!tru$wall.rhs] = tru$brk.rhs[!tru$wall.rhs] - tru$position[!tru$wall.rhs]
	
	
	truth[[tag]] = tru
}
truth = rbindlist(truth)


save(truth, file = "./truth/ibd_mrca.truth.RData")



