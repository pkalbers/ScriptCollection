#
# match results to truth
#

library(data.table)
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)

args = commandArgs(T)

prefix = args[1]             # file prefix
truth.file = args[2]             # truth data

load(sprintf("%s.RData", prefix))

load(truth.file)
load(sprintf("result.pairs.%s.RData", prefix))



truth = truth
del = which(truth$wall.lhs | truth$wall.rhs)
if (length(del) > 0) {
	truth = truth[-del, ]
}

truth.key = sprintf("%d %d %d", truth$index, truth$h0, truth$h1)


setup = as.data.frame(data)

del = which(setup$g.wall.lhs | setup$g.wall.rhs); if (length(del) > 0) setup = setup[-del, ]
del = which(setup$h.wall.lhs | setup$h.wall.rhs); if (length(del) > 0) setup = setup[-del, ]
del = which(setup$p.wall.lhs | setup$p.wall.rhs); if (length(del) > 0) setup = setup[-del, ]

setup$g.wall.lhs = NULL;  setup$g.wall.rhs = NULL
setup$h.wall.lhs = NULL;  setup$h.wall.rhs = NULL
setup$p.wall.lhs = NULL;  setup$p.wall.rhs = NULL

cat(" ", nrow(setup), "pairs ")

if (nrow(setup) == 0) next

setup.key = sprintf("%d %d %d", setup$index, setup$h0, setup$h1)

i = which(setup.key %in% truth.key)
j = which(truth.key %in% setup.key)

if (! identical(setup$index[i], truth$index[j])) stop("Mismatch: index")
if (! identical(setup$h0[i],    truth$h0[j])) stop("Mismatch: H0")
if (! identical(setup$h1[i],    truth$h1[j])) stop("Mismatch: H1")

sub = setup[i, ]
tru = truth[j, ]

cat(sprintf("--> matching %.1f%%\n", nrow(sub) / nrow(setup) * 100))

sub$t.lhs = tru$lhs
sub$t.rhs = tru$rhs

sub$t.dist.lhs = tru$dist.lhs
sub$t.dist.rhs = tru$dist.rhs

sub$t.brk.lhs = tru$brk.lhs
sub$t.brk.rhs = tru$brk.rhs

sub$t.brk.dist.lhs = tru$brk.dist.lhs
sub$t.brk.dist.rhs = tru$brk.dist.rhs

match = as.data.table(sub)



save(match, file = sprintf("result.match.%s.RData", prefix))


