
library(data.table)


tpl.generate = "PYTHON=/usr/local/var/homebrew/linked/python3/bin/python3.6;  ${PYTHON} ../pair2ms.py ../../vanilla.hdf5 %s %d %d"
#tpl.generate = "/users/mccarthy/pkalbers/miniconda2/bin/python ./pair2ms.py ../vanilla.hdf5 %s %d %d"
tpl.decode.m = "../decode_mod -m 0.0002 -r 0.0002 -p %d %s"


args = commandArgs(T)

file = args[1]

ccf = fread(file, header = T, stringsAsFactors = F)

M = fread("../../vanilla.marker.txt", header = T, stringsAsFactors = F)

pairs = data.table(mid = ccf$MarkerID, 
									 pos = M$Position[ ccf$MarkerID + 1 ],
									 id0 = (ccf$SampleID0 * 2) + ccf$Chr0,
									 id1 = (ccf$SampleID1 * 2) + ccf$Chr1,
									 shr = ccf$Shared)

pairs = unique(pairs)


psmc = list()

for (i in 1:nrow(pairs)) {
	
	mid = pairs$mid[i]
	pos = pairs$pos[i]
	id0 = pairs$id0[i]
	id1 = pairs$id1[i]
	shr = pairs$shr[i]
	
	tmp = tempfile()
	
	cat(mid, pos, id0, id1, shr, "...")
	
	system(sprintf(tpl.generate, tmp, id0, id1)) # generate psmc input data file
	
	if (is.na(file.size(tmp))) {
		cat(" ERROR(1)\n")
		next
	}
	
	out = system(sprintf(tpl.decode.m, pos, tmp), intern = T, ignore.stderr = T) # decode TMRCA using PSMC
	
	if (length(out) != 3) {
		cat(" ERROR(2)\n")
		next
	}
	
	unlink(tmp)
	
	tmp = list()
	tmp$times  = as.numeric(strsplit(out[1], split = " ", fixed = T)[[1]][-1])
	tmp$posterior = as.numeric(strsplit(out[3], split = " ", fixed = T)[[1]])
	
	target = as.numeric(strsplit(out[2], split = " ", fixed = T)[[1]][-1])
	
	if (pos != target) {
		cat(" ERROR(3)\n")
		next
	}
	
	tag = sprintf("%d %d %d %d %d", mid, pos, id0, id1, shr)
	psmc[[tag]] = tmp
	
	cat(" OK\n")
}



times = psmc[[1]]$times
check = sapply(psmc, function(x, y) {
	identical(x$times, y)
}, times)

if (all(check)) {
	psmc = lapply(psmc, function(x) {
		x$posterior
	})
}



save(times, psmc, file = sprintf("psmc_%s.RData", basename(file)))



