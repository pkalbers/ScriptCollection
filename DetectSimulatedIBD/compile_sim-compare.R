#
# match error results to originals
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

err.file = args[1] # eg. "result.pairs.history.generror_1000g.RData"
sim.file = args[2] # eg. "../result.pairs.history.RData"
out.file = args[3] # eg. "result.sim.history.generror_1000g.RData"


load(sim.file)

sim = data
rm(data)

load(err.file)

err = data
rm(data)


del = which(sim$g.wall.lhs | sim$g.wall.rhs); if (length(del) > 0) sim = sim[-del, ]
del = which(sim$h.wall.lhs | sim$h.wall.rhs); if (length(del) > 0) sim = sim[-del, ]
del = which(sim$p.wall.lhs | sim$p.wall.rhs); if (length(del) > 0) sim = sim[-del, ]

sim = as.data.frame(sim)
sim$g.wall.lhs = NULL;  sim$h.wall.lhs = NULL;  sim$p.wall.lhs = NULL;  
sim$g.wall.rhs = NULL;  sim$h.wall.rhs = NULL;  sim$p.wall.rhs = NULL;  

del = which(err$g.wall.lhs | err$g.wall.rhs); if (length(del) > 0) err = err[-del, ]
del = which(err$h.wall.lhs | err$h.wall.rhs); if (length(del) > 0) err = err[-del, ]
del = which(err$p.wall.lhs | err$p.wall.rhs); if (length(del) > 0) err = err[-del, ]

err = as.data.frame(err)
err$g.wall.lhs = NULL;  err$h.wall.lhs = NULL;  err$p.wall.lhs = NULL;  
err$g.wall.rhs = NULL;  err$h.wall.rhs = NULL;  err$p.wall.rhs = NULL;  


sim.key = sprintf("%d %d %d", sim$index, sim$g0, sim$g1)
err.key = sprintf("%d %d %d", err$index, err$g0, err$g1)


# intersection

key = intersect(sim.key, err.key)

sim.idx = match(key, sim.key)
err.idx = match(key, err.key)

inter.sim = sim[sim.idx, ]
inter.err = err[err.idx, ]

if (!identical(inter.sim$index, inter.err$index)) stop("index")
if (!identical(inter.sim$g0, inter.err$g0)) stop("g0")
if (!identical(inter.sim$g1, inter.err$g1)) stop("g1")


# diff sets

key = setdiff(sim.key, err.key)
sim.idx = match(key, sim.key)
diff.sim = sim[sim.idx, ]

key = setdiff(err.key, sim.key)
err.idx = match(key, err.key)
diff.err = err[err.idx, ]



save(inter.sim, diff.sim, inter.err, diff.err, file = out.file)


