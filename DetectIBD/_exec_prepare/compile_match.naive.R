#
# match results to truth
#

library(data.table)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20.GenErr_1000G" #args[1]
data.file = sprintf("./result.pairs.%s.RData", prefix)
true.file = "../result.truth.local.RData" #args[2]

load(data.file)
load(true.file)


data = as.data.frame(data)
truth = as.data.frame(truth)

data.key = sprintf("% 9d % 4d % 4d",  data$index,  data$g0,  data$g1)
true.key = sprintf("% 9d % 4d % 4d", truth$index, truth$g0, truth$g1)

data = data[order(data.key), ]
data.key = sort(data.key)

truth = truth[order(true.key), ]
true.key = sort(true.key)

key = intersect(data.key, true.key)

data.i = which(data.key %in% key)
true.i = which(true.key %in% key)

if (any(range(abs(data$index[data.i] - truth$index[true.i])) != 0)) stop("Mismatch: index")
if (any(range(abs(data$g0[data.i] - truth$g0[true.i])) != 0)) stop("Mismatch: G0")
if (any(range(abs(data$g1[data.i] - truth$g1[true.i])) != 0)) stop("Mismatch: G1")


dat = data[data.i, ]
tru = truth[true.i, ]

dat$true.time = tru$time
dat$true.lhs.pos = tru$lhs.position
dat$true.rhs.pos = tru$rhs.position
dat$true.lhs.idx = tru$lhs.index
dat$true.rhs.idx = tru$rhs.index
dat$true.wall = tru$wall


match = as.data.table(dat)



save(match, file = sprintf("result.match.%s.RData", prefix))


