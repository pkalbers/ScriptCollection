#
# match results to truth
#

library(data.table)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]
data.file = sprintf("./result.pairs_hmm_trunc.%s.RData", prefix)
true.file = sprintf("./result.match.%s.RData", prefix) #args[2]

load(data.file)
load(true.file)


data = as.data.frame(data)
truth = as.data.frame(match) #####

data.key = sprintf("% 9d % 4d % 4d",  data$index,  data$g0,  data$g1)
true.key = sprintf("% 9d % 4d % 4d", truth$index, truth$g0, truth$g1)

data = data[order(data.key), ]
data.key = sprintf("% 9d % 4d % 4d",  data$index,  data$g0,  data$g1)

truth = truth[order(true.key), ]
true.key = sprintf("% 9d % 4d % 4d", truth$index, truth$g0, truth$g1)

key = intersect(data.key, true.key)

data.i = which(data.key %in% key)
true.i = which(true.key %in% key)

if (any(range(abs(data$index[data.i] - truth$index[true.i])) != 0)) stop("Mismatch: index")
if (any(range(abs(data$g0[data.i] - truth$g0[true.i])) != 0)) stop("Mismatch: G0")
if (any(range(abs(data$g1[data.i] - truth$g1[true.i])) != 0)) stop("Mismatch: G1")


dat = data[data.i, ]
tru = truth[true.i, ]


tru$hmm.lhs = dat$hmm.lhs
tru$hmm.rhs = dat$hmm.rhs

tru$hmm.t8.lhs = dat$hmm.t8.lhs
tru$hmm.t8.rhs = dat$hmm.t8.rhs

tru$hmm.t32.lhs = dat$hmm.t32.lhs
tru$hmm.t32.rhs = dat$hmm.t32.rhs

tru$hmm.t128.lhs = dat$hmm.t128.lhs
tru$hmm.t128.rhs = dat$hmm.t128.rhs

tru$hmm.wall = dat$hmm.wall

match = as.data.table(tru)



save(match, file = sprintf("result.match_hmm_trunc.%s.RData", prefix))


