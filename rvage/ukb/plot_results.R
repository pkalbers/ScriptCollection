

library(data.table)
library(ggplot2)
library(ggthemes)


files = dir(pattern = ".+\\.age\\.sites\\.HMM\\.Rec\\.txt$")
files = dir(pattern = ".+\\.age\\.sites\\.HMM\\.Rec\\.HardBreaks\\.txt$")

d = lapply(files, function(file) {
	x = read.table(file, header = T, stringsAsFactors = F)
	if (nrow(x) == 0) return(NULL)
	x$Chr = as.numeric(sub(".+chr([0-9]+).+", "\\1", file))
	x$Impact = "-"
	if (grepl("neutral", file))
		x$Impact = "not annotated"
	if (grepl("impact", file))
		x$Impact = sub(".*impact_chr[0-9]+\\.([A-Z]+)\\..+", "\\1", file)
	if (grepl("pathog", file))
		x$Impact = "PATHOGENIC"
	if (grepl("pubmed", file))
		x$Impact = "PUBMED"
	if (grepl("random", file))
		x$Impact = "RANDOM"
	x
})
d = rbindlist(d)

#d = d[which(d$Fk >= 10),]

del = which(duplicated(sprintf("%d %d", d$Chr, d$MarkerID)))
if (length(del) > 0)
	d = d[-del,]


d$frq = d$Fk / (152249*2)
d$Frequency = d$Fk / (152249*2)

d = d[-(which(d$Impact == "-")),]




d = split(d, d$Chr)

RES = NULL
VEP = NULL

for (tag in names(d)) {
	cat(tag, "\n")
	
	sub = d[[tag]]
	
	load(sprintf("~/Research/UKbiobank/vep/chr%s.vcf_vep.RData", tag))
	
	mrk = fread(sprintf("../ukb_chr%s_c.marker.txt", tag))

	x = match(sub$MarkerID, mrk$MarkerID)
	mrk = mrk[x, ]
	
	
	sub$MarkerID = sprintf("%s:%d", tag, sub$MarkerID)
	mrk$MarkerID = sprintf("%s:%d", tag, mrk$MarkerID)
	
	
	mrk$RSID = sub(".+(rs.+)$", "\\1", mrk$Label)
	
	x = grep("rs", mrk$RSID, invert = T)
	if (length(x) > 0) {
		mrk$RSID[x] = "."
	}
	

	
	if (!identical(sub$MarkerID, mrk$MarkerID)) stop("???")
	
	
	x = which(vep$Uploaded_variation %in% mrk$RSID)
	vep = vep[x, ]
	
	
	x = match(vep$Uploaded_variation, mrk$RSID)
	vep = cbind(MarkerID = mrk$MarkerID[x], vep)
	
	
	sub = sub[order(sub$MarkerID), ]
	mrk = mrk[order(mrk$MarkerID), ]
	vep = vep[order(vep$MarkerID), ]
	
	if (!identical(sub$MarkerID, mrk$MarkerID)) stop("???")
	
	
	sub = as.data.frame(sub)
	mrk = as.data.frame(mrk)
	vep = as.data.frame(vep)
	
	sub = sub[, c("MarkerID","Fk","N_Shared","N_Others","PostMean","PostMode","PostMedian","PostCI025","PostCI975","Robust","Lower","Upper","Impact")]
	mrk = mrk[, c("Chromosome","Position","RSID","Alleles","AlleleCount0","AlleleCount1","AlleleCountX","GenotypeCount0","GenotypeCount1","GenotypeCount2","GenotypeCountX","RecombRate","GenDist" )]
	vep = vep[, c("MarkerID","Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extras")]
	
	res = cbind(sub, mrk)
	
	RES = rbind(RES, res)
	VEP = rbind(VEP, vep)
	
}


h.age = RES
h.vep = VEP

save(h.age, h.vep, file="../ukbb_data_vep.hardbreak.RData")






p = d[which(d$Impact == "PATHOGENIC" | d$Impact == "LOW" | d$Impact == "HIGH" | d$Impact == "MODERATE" | d$Impact == "MODIFIER"),]

p = p[which(p$Fk >= 75),]
p = p[which(p$Fk <= 155),]



p = rbindlist(lapply(split(p, p$Impact), function(x, n) x[sort(sample(1:nrow(x), n)), ], min(table(p$Impact))))




ggplot(p) +
	geom_point(aes(x = frq * 100, y = PostMean * 2 * 10000, colour = Impact), alpha = 0.75) +
	geom_smooth(aes(x = frq * 100, y = PostMean * 2 * 10000, colour = Impact), method = lm) +
	scale_x_log10() +
	scale_y_log10(breaks=c(10^(0:5)), labels=sprintf("%d", 10^(0:5))) +
	scale_color_brewer(palette = "Set1") +
	theme_few() +
	ylab("Estimated age") + xlab("Allele frequency (%)")



rm = split(p, p$Impact)
rm = lapply(rm, function(x) { 
	x = x[order(x$Fk), ]
	n = 50 # round(nrow(x) / 10)
	data.table(frq = rollmean(x$frq, n), res = rollmean(x$PostMean, n), Impact = x$Impact[1]) 
} )
rm = rbindlist(rm)


ggplot(rm) +
	geom_point(aes(x = frq * 100, y = res * 2 * 10000, colour = Impact), alpha = 0.75) +
	#geom_smooth(aes(x = fk, y = res * 2 * 10000, colour = imp), method = lm) +
	scale_x_log10() +
	scale_y_log10(breaks=c(10^(0:5)), labels=sprintf("%d", 10^(0:5))) +
	scale_color_brewer(palette = "Set1") +
	theme_few() +
	ylab("Estimated age") + xlab("Allele frequency (%)")



ggplot(p) +
	geom_freqpoly(aes(frq, colour = Impact)) +
	scale_y_log10()


