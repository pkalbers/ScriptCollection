#
# compile all results
#

library(data.table)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20.GenErr_1000G" #args[1]

load(sprintf("data.%s.RData", prefix))



wall.lhs = 50000
wall.rhs = POS[length(POS)] - 50000


tmp.files = dir(pattern = "^tmp_hmm_trunc\\.pairs\\..+\\.([0-9]+)\\.RData$", path = "./pairs_hmm", full.names = T) #####
#res.files = dir(pattern = "^result.+\\.pairs\\..+\\.([0-9]+)\\.RData$", path = "./pairs_hmm", full.names = T) #####

data = list()

for (tmp.file in tmp.files) {
	pack = sub(".+\\.([0-9]+)\\.RData", "\\1", tmp.file)
	
	cat(pack)
	
	if (file.exists(sub("tmp", "result", tmp.file)))
		load(sub("tmp", "result", tmp.file))
	else
		load(tmp.file)
	
	# removing failed detections
	
	fail = which(is.na(pair$hmm.lhs) | is.na(pair$hmm.rhs) | (pair$hmm.lhs == 0 & pair$hmm.rhs == 0))
	if (length(fail) > 0) {
		cat(sprintf(" [Failed:  %.3f%%]", length(fail) / nrow(pair) * 100))
		pair = pair[-fail, ]
		if (nrow(pair) == 0) {
			next
		}
	}
	
	
	# reaches end of chromosome
	
	pair$hmm.wall = F;  wall = which(pair$hmm.lhs < wall.lhs | pair$hmm.rhs > wall.rhs);  if (length(wall) > 0) pair$hmm.wall[wall] = T
	
	
	cat("\n")
	
	data[[pack]] = pair
}
data = rbindlist(data)


save(data, file = sprintf("result.pairs_hmm_trunc.%s.RData", prefix))



