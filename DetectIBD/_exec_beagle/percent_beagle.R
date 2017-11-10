#
# plot detected vs truth
#

library(data.table)
library(ggplot2)
library(ggthemes)
library(colorspace)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1]

load(sprintf("data.%s.RData", prefix))

DIST = c(0, cumsum(RATE[-(length(RATE))] * diff(POS) * 1e-6)) ### genetic distance
names(DIST) = as.character(POS)


load("./result.truth.local.RData")
load(sprintf("result.match.%s.RData", prefix))
load(sprintf("result.beagle.%s.RData", prefix))


(nrow(match) / nrow(truth)) * 100
(nrow(match.H) / nrow(truth)) * 100
(nrow(match.P) / nrow(truth)) * 100


key = sprintf("%d %d %d %d %d", truth$fk, truth$g0, truth$g1, truth$lhs.index, truth$rhs.index)
key = unique(key)


match$key = sprintf("%d %d %d %d %d", match$fk, match$g0, match$g1, match$, match$true.rhs.idx)
# match = match[which(match$key %in% key), ]


match.H$key = sprintf("%d %d %d %d %d", match.H$fk, match.H$g0, match.H$g1, match.H$beagle.lhs.index, match.H$beagle.rhs.index)
match.P$key = sprintf("%d %d %d %d %d", match.P$fk, match.P$g0, match.P$g1, match.P$beagle.lhs.index, match.P$beagle.rhs.index)

# match.H = match.H[which(match.H$key %in% key), ]
# match.P = match.P[which(match.P$key %in% key), ]


m = match
del = which(duplicated(match$key))
if (length(del) > 0) {
	m = match[-del, ]
}

h = match.H
del = which(duplicated(match.H$key))
if (length(del) > 0) {
	h = match.H[-del, ]
}

p = match.P
del = which(duplicated(match.P$key))
if (length(del) > 0) {
	p = match.P[-del, ]
}




nrow(m)
nrow(h)
nrow(p)

(nrow(m) / length(key)) * 100
(nrow(h) / length(key)) * 100
(nrow(p) / length(key)) * 100

(table(h$fk) / table(m$fk)) * 100
(table(p$fk) / table(m$fk)) * 100













