
library(data.table)


tpl.generate = "PYTHON=/usr/local/var/homebrew/linked/python3/bin/python3.6;  ${PYTHON} ../hap2psmc.py %s %s"
#tpl.generate = "/users/mccarthy/pkalbers/miniconda2/bin/python ../hap2psmc.py %s %s"
tpl.decode.m = "../decode_mod -m 0.00024 -r 0.00034 -p %d %s"  # 1.2e-8 * 20000 , 1.7e-8 * 20000 ( (108.266934 / 100) / (62949445 - 61795) )


args = commandArgs(T)

file = args[1]

ccf = fread(file, header = T, stringsAsFactors = F)

load("../snp_1000G_chr20_rle_data.RData")
POS=pos
TAG=tag

pairs = data.table(mid = ccf$MarkerID, 
									 pos = POS[ ccf$MarkerID + 1 ],
									 id0 = (ccf$SampleID0 * 2) + ccf$Chr0,
									 id1 = (ccf$SampleID1 * 2) + ccf$Chr1,
									 shr = ccf$Shared)

pairs = unique(pairs)


pos.str = paste(POS, collapse = ",")


psmc = list()

for (i in 1:nrow(pairs)) {
	
	mid = pairs$mid[i]
	pos = pairs$pos[i]
	id0 = pairs$id0[i]
	id1 = pairs$id1[i]
	shr = pairs$shr[i]
	
	
	cat(mid, pos, id0, id1, shr, "...")
	
	
	# get haplotypes
	
	v0 = inverse.rle(H[[ TAG[ id0+1 ] ]])
	v1 = inverse.rle(H[[ TAG[ id1+1 ] ]])
	
	if (length(v0) != length(POS) || length(v1) != length(POS)) {
		cat(" DATA ERROR (1)")
		next
	}
	
	if (v0[mid+1] + v1[mid+1] != shr + 1) {
		cat(" DATA ERROR (2)")
		next
	}
	
	s0 = paste(v0, collapse = '')
	s1 = paste(v1, collapse = '')
	
	tmp1 = tempfile()
	tmp2 = tempfile()
	
	cat(c(pos.str, s0, s1), sep = "\n", file = tmp1)
	
	
	system(sprintf(tpl.generate, tmp1, tmp2)) # generate psmc input data file
	
	if (is.na(file.size(tmp2))) {
		cat(" ERROR(1)\n")
		next
	}
	
	out = system(sprintf(tpl.decode.m, pos, tmp2), intern = T, ignore.stderr = T) # decode TMRCA using PSMC
	
	if (length(out) != 3) {
		cat(" ERROR(2)\n")
		next
	}
	
	unlink(tmp1)
	unlink(tmp2)
	
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



