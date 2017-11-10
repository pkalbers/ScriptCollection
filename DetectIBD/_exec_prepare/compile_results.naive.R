#
# compile all results
#

library(data.table)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1] 

load(sprintf("data.%s.RData", prefix))



wall.lhs = 50000
wall.rhs = POS[length(POS)] - 50000


res.files = dir(pattern = "^result.+\\.pairs\\..+\\.([0-9]+)\\.RData$", path = "./pairs", full.names = T)
res.packs = unique(sub("^.+\\.([0-9]+)\\.RData$", "\\1", res.files))

data = list()

for (pack in res.packs) {
	cat(pack)
	
	res.g = sprintf("./pairs/result_G.pairs.%s.%s.RData", prefix, pack)
	res.h = sprintf("./pairs/result_H.pairs.%s.%s.RData", prefix, pack)
	res.p = sprintf("./pairs/result_P.pairs.%s.%s.RData", prefix, pack)
	
	if (!file.exists(res.g) || !file.exists(res.h) || !file.exists(res.p)) {
		next
	}
	
	load(res.g)
	tmp = pair
	load(res.h)
	if (!identical(tmp$index, pair$index)) stop("Mismatch!!!")
	tmp = cbind(tmp, h.lhs = pair$h.lhs, h.rhs = pair$h.rhs)
	load(res.p)
	if (!identical(tmp$index, pair$index)) stop("Mismatch!!!")
	tmp = cbind(tmp, p.lhs = pair$p.lhs, p.rhs = pair$p.rhs)
	
	pair = tmp

	# removing failed detections
	
	fail.g = which(is.na(pair$g.lhs) | is.na(pair$g.rhs) | (pair$g.lhs == 0 & pair$g.rhs == 0))
	fail.h = which(is.na(pair$h.lhs) | is.na(pair$h.rhs) | (pair$h.lhs == 0 & pair$h.rhs == 0))
	fail.p = which(is.na(pair$p.lhs) | is.na(pair$p.rhs) | (pair$p.lhs == 0 & pair$p.rhs == 0))
	fail = unique(c(fail.g, fail.h, fail.p))
	if (length(fail) > 0) {
		cat(sprintf(" [Failed:  %.3f%%]", length(fail) / nrow(pair) * 100))
		pair = pair[-fail, ]
		if (nrow(pair) == 0) {
			next
		}
	}
	
	
	# reaches end of chromosome
	
	pair$g.wall = F;  wall = which(pair$g.lhs < wall.lhs | pair$g.rhs > wall.rhs);  if (length(wall) > 0) pair$g.wall[wall] = T
	pair$h.wall = F;  wall = which(pair$h.lhs < wall.lhs | pair$h.rhs > wall.rhs);  if (length(wall) > 0) pair$h.wall[wall] = T
	pair$p.wall = F;  wall = which(pair$p.lhs < wall.lhs | pair$p.rhs > wall.rhs);  if (length(wall) > 0) pair$p.wall[wall] = T
	

	cat("\n")
	
	data[[pack]] = pair
}
data = rbindlist(data)


save(data, file = sprintf("result.pairs.%s.RData", prefix))



