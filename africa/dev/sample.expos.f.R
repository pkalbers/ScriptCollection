

library(data.table)


mt = fread("~/Research/rvage/africa/truH.marker.txt")
me = fread("~/Research/rvage/africa/errH.marker.txt")


i = which(mt$AlleleCount1 == mt$GenotypeCount1)

k0 = intersect(which(mt$AlleleCount1 > 1), which(me$AlleleCount1 > 1))
k1 = intersect(which(mt$AlleleCount1 <= 61), which(me$AlleleCount1 <= 61))
k = intersect(k0, k1)

k = intersect(k, i)


f.id = intersect(which(mt$AlleleCount1 == me$AlleleCount1), k)
f.fn = intersect(which(mt$AlleleCount1 >  me$AlleleCount1), k)
f.fp = intersect(which(mt$AlleleCount1 <  me$AlleleCount1), k)


tmp = split(f.id, me$AlleleCount1[f.id])
tmp = split(f.fn, me$AlleleCount1[f.fn])
tmp = split(f.fp, me$AlleleCount1[f.fp])


tmp = lapply(tmp, function(x) {
	if (length(x) < 100) return(x)
	sort(sample(x, 100))
})
tmp = unlist(tmp)

del = which(me$AlleleCount1[tmp] > 51)
if (length(del) > 0) tmp = tmp[-del]

table(me$AlleleCount1[tmp])


p.id = me$Position[tmp]
p.fn = me$Position[tmp]
p.fp = me$Position[tmp]


cat(p.id, sep = "\n", file = "~/Research/rvage/africa/expos.f_id.txt")
cat(p.fn, sep = "\n", file = "~/Research/rvage/africa/expos.f_fn.txt")
cat(p.fp, sep = "\n", file = "~/Research/rvage/africa/expos.f_fp.txt")
