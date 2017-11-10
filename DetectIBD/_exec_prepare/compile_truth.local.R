#
# compile truth
#

library(data.table)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("data.%s.RData", prefix))


wall.lhs = 50000
wall.rhs = POS[length(POS)] - 50000

tru.files = dir(pattern = "^local_ibd\\.pairs\\..+\\.txt$", path = "./_truth", full.names = T)


truth = list()

for (tru.file in tru.files) {
	tru = read.table(tru.file, header = T)
	tru = as.data.table(tru)
	
	tag = sub("^.+\\.([0-9]+)\\..+$", "\\1", tru.file)
	
	cat(tag, "\n")
	
	# reaches end of chromosome
	tru$wall = F; wall = which(tru$lhs.position < wall.lhs | tru$rhs.position > wall.rhs); if (length(wall) > 0) tru$wall[wall] = T
	
	tru$lhs.index = tru$lhs.index + 1
	tru$rhs.index = tru$rhs.index + 1
	
	truth[[tag]] = tru
}
truth = rbindlist(truth)


save(truth, file = "./result.truth.local.RData")



