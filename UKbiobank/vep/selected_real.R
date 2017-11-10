
library(data.table)

real = NULL

for (i in 1:22) {
	cat(i, "\n")
	
	load(sprintf("pruned.chr%d.pub_pat.RData", i))
	
	k = which(pat$impact == "HIGH")
	
	if (length(k) > 0) {
		pat = pat[k, ]
		
		k = grep("^[0-9,]+$", pat$pubmed)
		
		if (length(k) > 0) {
			pat = pat[k, ]
			
			k = grep(".+rs.+", pat$Label)
			
			if (length(k) > 0) {
				pat = pat[k, ]
				
				del = which(duplicated(pat$MarkerID))
				if (length(del) > 0)
					pat = pat[-del, ]
				
				real = rbind(real, pat)
				
			}
		}
	}
}


save(real, file = "selected.real.RData")


stop()


load("selected.real.RData")


out = data.table(ref = real$MarkerID, chr = real$Chromosome, pos = real$Position)

out = rbind(out, data.table(ref = "HCM", chr = 1 , pos = 236911044))

out = as.data.frame(out)


write.table(out, file = "selected_real.txt", append = F, sep = " ", row.names = F, col.names = F)

