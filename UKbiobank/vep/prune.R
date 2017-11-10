

library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	load(sprintf("chr%d.vcf_vep.RData", i))
	
	m = fread(sprintf("ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	del = which(m$AlleleCount1 != m$GenotypeCount1)
	if (length(del) > 0)
		m = m[-del, ]
	
	del = which(m$AlleleCount1 < 2)
	if (length(del) > 0)
		m = m[-del, ]
	
	size = m$AlleleCount0[1] + m$AlleleCount1[1] + m$AlleleCountX[1]
	
	del = which(m$AlleleCount1 > size*0.005) # FK LIMIT !!!
	if (length(del) > 0)
		m = m[-del, ]
	
	m$rsid = sub(".+_(rs[0-9]+)$", "\\1", m$Label)
	
	
	#net = which(!m$rsid %in% unique(vep$Uploaded_variation))
	#net = m[net, ]
	#save(net, file = sprintf("pruned.chr%d.net.RData", i))
	#next
	
	
	del = grep("^rs.+", m$rsid, invert = T)
	if (length(del) > 0)
		m = m[-del, ]
	
	x = match(m$Position, vcf$pos)
	
	del = which(is.na(x))
	if (length(del) > 0) {
		m = m[-del, ]
		x = x[-del]
	}
	
	
	m$vcf_pos = vcf$pos[x]
	m$vcf_pos38 = vcf$POS38[x]
	m$vcf_alleles = sprintf("%s,%s", vcf$a0[x], vcf$a1[x])
	
	del = which(m$vcf_pos != m$Position)
	if (length(del) > 0)
		m = m[-del, ]
	
	del = which(m$vcf_alleles != m$Alleles)
	if (length(del) > 0)
		m = m[-del, ]
	
	
	x = match(vep$Uploaded_variation, m$rsid)
	
	vep = cbind(vep, m[x, ])
	
	del = which(is.na(vep$MarkerID))
	if (length(del) > 0)
		vep = vep[-del, ]
	
	
	vep$impact = sub("^IMPACT=([^;]+);.+", "\\1", vep$Extras)
	
	
	pub = grep("PUBMED", vep$Extras)
	if (length(pub) > 0) {
		pub = vep[pub, ]
		pub$pubmed = sub(".*PUBMED=([0-9,]+);?.*", "\\1", pub$Extras)
	} else {
		pub = NULL
	}
	
	pat = grep("patho", vep$Extras)
	if (length(pat) > 0) {
		pat = vep[pat, ]
		pat$pubmed = sub(".*PUBMED=([0-9,]+);?.*", "\\1", pat$Extras)
	} else {
		pat = NULL
	}
	
	save(pub, pat, file = sprintf("pruned.chr%d.pub_pat.RData", i))
	
	
	#save(vep, file = sprintf("pruned.chr%d.vep.RData", i))
}

















