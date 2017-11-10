

library(data.table)


M = list()

for (i in 1:22) {
	
	m = fread(sprintf("../ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	m$Chr = i
	
	M[[as.character(i)]] = m
}



A = fread("UKBB_HGMD_vars_ann.txt", header = T, stringsAsFactors = F)

del = which(is.na(A$UKBBIDBGEN))
if (length(del) > 0)
	A = A[-del, ]


A = split(A, A$CHR)


out = NULL

for (tag in names(A)) {
	print(tag)
	
	m = M[[tag]]
	a = A[[tag]]
	
	del = which(m$AlleleCount1 < 2)
	if (length(del) > 0)
		m = m[-del, ]
	
	x = match(a$POS_HG19, m$Position)
	
	if (length(x) > 0) {
		m = m[x, ]
		
		if (nrow(m) > 0) {
		
			out = rbind(out, cbind(m, a))
				
		}
	}
	
}


del = which(is.na(out$MarkerID))
if (length(del) > 0)
	out = out[-del, ]


del = which(out$AlleleCount1 > 1515)
if (length(del) > 0)
	out = out[-del, ]


hgmd = out


save(hgmd, file = "hgmd_ukbb.RData")







z = split(hgmd, hgmd$Chromosome)

for (tag in names(z)) {
	print(tag)
	
	x = z[[tag]]
	x = sprintf("%d|%d", x$MarkerID, x$Position)
	cat(x, sep = "\n", file = sprintf("hgmd_chr%s.txt", tag))
}





###
stop()


B = fread("HGMD_not_in_18kexome_230513_GMcV_CancerAnnot_manifest_ann.txt", header = T, stringsAsFactors = F)


# del = which(is.na(B$ukbb.rsid))
# if (length(del) > 0)
# 	B = B[-del, ]

#B$CHR = as.numeric(sub("^([0-9]+):.+$", "\\1", B$ukbb.snpid))

del = which(is.na(B$CHR))
if (length(del) > 0)
	B = B[-del, ]


B = split(B, B$CHR)


res = NULL

for (tag in names(B)) {
	print(tag)

	m = M[[tag]]
	b = B[[tag]]
	
	del = which(m$AlleleCount1 < 2)
	if (length(del) > 0)
		m = m[-del, ]
	
	x = match(b$UKB_WCSG_POS, m$Position)
	
	if (length(x) > 0) {
		m = m[x, ]
		
		if (nrow(m) > 0) {
			
			res = rbind(res, cbind(m, b))
			
		}
		
	}
}


del = which(is.na(res$MarkerID))
if (length(del) > 0)
	res = res[-del, ]


del = which(res$AlleleCount1 > 1515)
if (length(del) > 0)
	res = res[-del, ]



length(unique(c(out$MarkerID, res$MarkerID)))









