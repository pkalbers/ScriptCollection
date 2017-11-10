#
# run global IBD identification
#

library(data.table)
library(rPython)


#args = commandArgs(T)

prefix = "OutOfAfricaHapMap20" #args[1] 

sample.n = 1000

hist.file = "../OutOfAfricaHapMap20.hdf5"           ###### path to simulation history file
load.file = "../_exec_msprime/msprime.fetch_ibd.py" ###### path to fetch file

load(sprintf("../data.%s.RData", prefix))



fetch.ibd = function(h0, h1, return.raw = F, history = hist.file, load = load.file) {
	python.load(load)
	ibd = python.call("MRCA", h0, h1, history)
	ibd = data.table(mrca = ibd$mrca, 
									 time = ibd$time, 
									 lhs.position = ibd$lhs.position, 
									 rhs.position = ibd$rhs.position, 
									 lhs.index = ibd$lhs.index, 
									 rhs.index = ibd$rhs.index)
	
	if (return.raw) {
		return (ibd)
	}
	
	above = which(ibd$mrca[-(nrow(ibd))] != ibd$mrca[-1])
	below = c(1, above + 1)
	above = c(above, nrow(ibd))
	
	data.table(mrca = ibd$mrca[below], 
						 time = ibd$time[below], 
						 lhs.position = ibd$lhs.position[below], 
						 rhs.position = ibd$rhs.position[above], 
						 lhs.index = ibd$lhs.index[below], 
						 rhs.index = ibd$rhs.index[above])
}



run.pair = function(g0, g1) {
	pair = sprintf("%d %d", g0, g1)
	
	cat(pair, " ")
	
	A0 = if (g0 %% 2 != 0) g0 * 2 else g0 * 2 - 1
	A1 = if (g0 %% 2 != 0) g0 * 2 - 1 else g0 * 2
	B0 = if (g1 %% 2 != 0) g1 * 2 else g1 * 2 - 1
	B1 = if (g1 %% 2 != 0) g1 * 2 - 1 else g1 * 2
	
	tmp = list()
	
	a0.a1 = sprintf("%d %d", A0, A1)
	tmp[[a0.a1]] = fetch.ibd(A0, A1)
	cat(".")
	
	a0.b0 = sprintf("%d %d", A0, B0)
	tmp[[a0.b0]] = fetch.ibd(A0, B0)
	cat(".")
	
	a0.b1 = sprintf("%d %d", A0, B1)
	tmp[[a0.b1]] = fetch.ibd(A0, B1)
	cat(".")
	
	a1.b0 = sprintf("%d %d", A1, B0)
	tmp[[a1.b0]] = fetch.ibd(A1, B0)
	cat(".")
	
	a1.b1 = sprintf("%d %d", A1, B1)
	tmp[[a1.b1]] = fetch.ibd(A1, B1)
	cat(".")
	
	b0.b1 = sprintf("%d %d", B0, B1)
	tmp[[b0.b1]] = fetch.ibd(B0, B1)
	cat(".\n")
	
	out = list()
	out[[pair]] = tmp
	out
}



ibd = c()

fki = FKI[sample(nrow(FKI)), ]

key = sprintf("%d %d", fki$g0, fki$g1)
del = which(duplicated(key))
fki = fki[-del, ]

save.file = sprintf("global_ibd.%s.RData", paste(sample(c(letters, LETTERS, as.character(0:9)), 16), collapse = ""))

for (i in 1:nrow(fki)) {
	g0 = fki$g0[i]
	g1 = fki$g1[i]
	
	tmp = run.pair(g0, g1)
	
	ibd = c(ibd, tmp)
	
	if (i == sample.n) {
		cat("Saving data\n")
		save(ibd, file = save.file)
		break
	}
}












