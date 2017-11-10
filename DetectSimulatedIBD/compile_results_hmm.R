#
# compile all results
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)


args = commandArgs(T)

prefix = args[1]             # file prefix

load(sprintf("%s.RData", prefix))


index = 1:length(POS)
index.names = sprintf("%.5f", POS)
names(index) = index.names


wall.lhs = POS[10]
wall.rhs = POS[length(POS)-10]


res.files = dir(pattern = "^result_hmm\\.pairs_hmm\\..+\\.([0-9]+)\\.RData$", path = "./pairs_hmm", full.names = T)

data = list()

for (res.file in res.files) {
	load(res.file)
	
	pack = sub("^.+\\.([0-9]+)\\..+$", "\\1", res.file)
	
# 	setup = sub("^([^\\/]*)\\/?pairs\\/.+", "\\1", res.file)
# 	if (grepl("^generror.+", setup)) setup = gsub("_", " ", sub("^generror_(.+)", "\\1", setup))
# 	if (setup == "") setup = "simulated"
# 	setup = toupper(setup)
# 	pair$setup = setup
	
	cat(pack)
	
	
	# removing failed detections
	
	fail = which(is.na(pair$m.lhs) | is.na(pair$m.rhs) | (pair$m.lhs == 0 & pair$m.rhs == 0))
	if (length(fail) > 0) {
		cat(sprintf(" [Failed:  %.3f%%]", length(fail) / nrow(pair) * 100))
		pair = pair[-fail, ]
		if (nrow(pair) == 0) {
			next
		}
	}
	
	i = which(pair$m.lhs == 0); if (length(i) > 0) pair$m.lhs[i] = pair$position[i]
	i = which(pair$m.rhs == 0); if (length(i) > 0) pair$m.rhs[i] = pair$position[i]
	
	
	# reaches end of chromosome
	
	pair$m.wall.lhs = F; wall = which(round(pair$m.lhs) < wall.lhs); if (length(wall) > 0) pair$m.wall.lhs[wall] = T
	pair$m.wall.rhs = F; wall = which(round(pair$m.rhs) > wall.rhs); if (length(wall) > 0) pair$m.wall.rhs[wall] = T
	
	
	# distance to focal position
	
	pair$m.dist.lhs = pair$m.lhs - pair$position
	pair$m.dist.rhs = pair$m.rhs - pair$position
	
	
	# breakpoint positions
	
	prebreak.m.lhs = sprintf("%.5f", pair$m.lhs[!pair$m.wall.lhs])
	prebreak.m.rhs = sprintf("%.5f", pair$m.rhs[!pair$m.wall.rhs])
	
	if (!all(prebreak.m.lhs %in% index.names) || !all(prebreak.m.rhs %in% index.names)) {
		stop("Cannot match all breakpoint positions")
	}
	
	pair$m.brk.lhs = NA;  pair$m.brk.lhs[!pair$m.wall.lhs] = POS[ index[prebreak.m.lhs] - 1 ]
	pair$m.brk.rhs = NA;  pair$m.brk.rhs[!pair$m.wall.rhs] = POS[ index[prebreak.m.rhs] + 1 ]
	
	
	# breakpoint distances
	
	pair$m.brk.dist.lhs = NA;  pair$m.brk.dist.lhs[!pair$m.wall.lhs] = pair$m.brk.lhs[!pair$m.wall.lhs] - pair$position[!pair$m.wall.lhs]
	pair$m.brk.dist.rhs = NA;  pair$m.brk.dist.rhs[!pair$m.wall.rhs] = pair$m.brk.rhs[!pair$m.wall.rhs] - pair$position[!pair$m.wall.rhs]
	
	
	cat("\n")
	
	data[[res.file]] = pair
}
data = rbindlist(data)



save(data, file = sprintf("result_hmm.pairs_hmm.%s.RData", prefix))



