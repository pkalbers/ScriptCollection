#
# compile all results
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)


prefix = "1000GP_Phase3_chr20"


# get positions
POS = read.table(sprintf("%s.legend", prefix), header = T, stringsAsFactors = F)
POS = POS$position



index = 1:length(POS)
index.names = sprintf("%.5f", POS)
names(index) = index.names


wall.lhs = POS[10]
wall.rhs = POS[length(POS)-10]


res.files = dir(pattern = "^result\\.pairs\\..+\\.([0-9]+)\\.RData$", path = "./pairs", full.names = T)

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
	
	fail.g = which(is.na(pair$g.lhs) | is.na(pair$g.rhs) | (pair$g.lhs == 0 & pair$g.rhs == 0))
	fail.h = which(is.na(pair$h.lhs) | is.na(pair$h.rhs) | (pair$h.lhs == 0 & pair$h.rhs == 0))
	fail = unique(c(fail.g, fail.h))
	if (length(fail) > 0) {
		cat(sprintf(" [Failed:  %.3f%%]", length(fail) / nrow(pair) * 100))
		pair = pair[-fail, ]
		if (nrow(pair) == 0) {
			next
		}
	}
	
	
	# reaches end of chromosome
	
	pair$g.wall.lhs = F; wall = which(round(pair$g.lhs) < wall.lhs); if (length(wall) > 0) pair$g.wall.lhs[wall] = T
	pair$g.wall.rhs = F; wall = which(round(pair$g.rhs) > wall.rhs); if (length(wall) > 0) pair$g.wall.rhs[wall] = T
	
	pair$h.wall.lhs = F; wall = which(round(pair$h.lhs) < wall.lhs); if (length(wall) > 0) pair$h.wall.lhs[wall] = T
	pair$h.wall.rhs = F; wall = which(round(pair$h.rhs) > wall.rhs); if (length(wall) > 0) pair$h.wall.rhs[wall] = T
	
	
	# distance to focal position
	
	pair$g.dist.lhs = pair$g.lhs - pair$position
	pair$g.dist.rhs = pair$g.rhs - pair$position
	
	pair$h.dist.lhs = pair$h.lhs - pair$position
	pair$h.dist.rhs = pair$h.rhs - pair$position
	
	
	# breakpoint positions
	
	prebreak.g.lhs = sprintf("%.5f", pair$g.lhs[!pair$g.wall.lhs])
	prebreak.g.rhs = sprintf("%.5f", pair$g.rhs[!pair$g.wall.rhs])
	
	if (!all(prebreak.g.lhs %in% index.names) || !all(prebreak.g.rhs %in% index.names)) {
		stop("Cannot match all breakpoint positions")
	}
	
	prebreak.h.lhs = sprintf("%.5f", pair$h.lhs[!pair$h.wall.lhs])
	prebreak.h.rhs = sprintf("%.5f", pair$h.rhs[!pair$h.wall.rhs])
	
	if (!all(prebreak.h.lhs %in% index.names) || !all(prebreak.h.rhs %in% index.names)) {
		stop("Cannot match all breakpoint positions")
	}
	
	pair$g.brk.lhs = NA;  pair$g.brk.lhs[!pair$g.wall.lhs] = POS[ index[prebreak.g.lhs] - 1 ]
	pair$g.brk.rhs = NA;  pair$g.brk.rhs[!pair$g.wall.rhs] = POS[ index[prebreak.g.rhs] + 1 ]
	
	pair$h.brk.lhs = NA;  pair$h.brk.lhs[!pair$h.wall.lhs] = POS[ index[prebreak.h.lhs] - 1 ]
	pair$h.brk.rhs = NA;  pair$h.brk.rhs[!pair$h.wall.rhs] = POS[ index[prebreak.h.rhs] + 1 ]

	
	
	# breakpoint distances
	
	pair$g.brk.dist.lhs = NA;  pair$g.brk.dist.lhs[!pair$g.wall.lhs] = pair$g.brk.lhs[!pair$g.wall.lhs] - pair$position[!pair$g.wall.lhs]
	pair$g.brk.dist.rhs = NA;  pair$g.brk.dist.rhs[!pair$g.wall.rhs] = pair$g.brk.rhs[!pair$g.wall.rhs] - pair$position[!pair$g.wall.rhs]
	
	pair$h.brk.dist.lhs = NA;  pair$h.brk.dist.lhs[!pair$h.wall.lhs] = pair$h.brk.lhs[!pair$h.wall.lhs] - pair$position[!pair$h.wall.lhs]
	pair$h.brk.dist.rhs = NA;  pair$h.brk.dist.rhs[!pair$h.wall.rhs] = pair$h.brk.rhs[!pair$h.wall.rhs] - pair$position[!pair$h.wall.rhs]
	
	
	cat("\n")
	
	data[[res.file]] = pair
}
data = rbindlist(data)


save(data, file = sprintf("result.pairs.%s.RData", prefix))



