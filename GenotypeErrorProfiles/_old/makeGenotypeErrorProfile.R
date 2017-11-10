#
# make GT error profile
#

args = commandArgs(T)

true.file = args[1]
call.file = args[2]
save.file = args[3]



cat("Loading data ... ")
load(true.file)
load(call.file)
cat("OK\n")


p$true.is.confident = NA


cat("Distributing data ... ")
match = intersect(rownames(p), c(names(truth), conf.hom.ref))
p.matched = p[match, ]
p.unmatch = p[-(which(rownames(p) %in% match)), ]
cat("OK\n")



cat("Matching true genotypes ... ")
match = intersect(rownames(p.matched), names(truth))
p.truth = p.matched[match, ]
p.truth$true.gt = truth[match]
p.truth$true.is.typed = T
p.truth$true.is.confident = F
cat("OK\n")



cat("Matching masked true genotypes ... ")
i = which(rownames(p.truth) %in% truth.mask)
p.truth$true.is.confident[i] = T
cat("OK\n")



cat("Matching assumed hom-ref genotypes ... ")
match = intersect(rownames(p.matched), conf.hom.ref)
p.cnfhr = p.matched[match, ]
p.cnfhr$true.gt = 0
p.cnfhr$true.is.typed = F
p.cnfhr$true.is.confident = T
cat("OK\n")



cat("Compiling ... ")
rm(p)
p = rbind(p.truth, p.cnfhr, p.unmatch)
cat("OK\n")


cat("Saving ... ")
save(p, file = save.file)
cat("OK\n")

