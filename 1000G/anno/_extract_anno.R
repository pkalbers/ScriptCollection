


VEP=c("transcript_ablation" = "high",
			"splice_acceptor_variant" = "high",
			"splice_donor_variant" = "high",
			"stop_gained" = "high",
			"frameshift_variant" = "high",
			"stop_lost" = "high",
			"start_lost" = "high",
			"transcript_amplification" = "high",
			"inframe_insertion" = "moderate",
			"inframe_deletion" = "moderate",
			"missense_variant" = "moderate",
			"protein_altering_variant" = "moderate",
			"splice_region_variant" = "low",
			"incomplete_terminal_codon_variant" = "low",
			"stop_retained_variant" = "low",
			"synonymous_variant" = "low",
			"coding_sequence_variant" = "modifier",
			"mature_miRNA_variant" = "modifier",
			"5_prime_UTR_variant" = "modifier",
			"3_prime_UTR_variant" = "modifier",
			"non_coding_transcript_exon_variant" = "modifier",
			"intron_variant" = "modifier",
			"NMD_transcript_variant" = "modifier",
			"non_coding_transcript_variant" = "modifier",
			"upstream_gene_variant" = "modifier",
			"downstream_gene_variant" = "modifier",
			"TFBS_ablation" = "modifier",
			"TFBS_amplification" = "modifier",
			"TF_binding_site_variant" = "modifier",
			"regulatory_region_ablation" = "moderate",
			"regulatory_region_amplification" = "modifier",
			"feature_elongation" = "modifier",
			"regulatory_region_variant" = "modifier",
			"feature_truncation" = "modifier",
			"intergenic_variant" = "modifier")


lookup = names(VEP)



library(data.table)

files = dir(pattern = "^rsids\\.chr([0-9]+)\\.txt$")

for (file in files) {

	print(file)

	CHR = sub("^rsids\\.chr([0-9]+)\\.txt$", "\\1", file)

anno = fread(file, header = F, stringsAsFactors = F)

names(anno) = c("chr", "pos", "rsid", "a0", "a1", "qual", "pass", "pred")

anno$SOterm = ""
anno$impact = ""

for (look in lookup) {
	print(look)

	imp = as.character(VEP[which(lookup == look)])

	x = which(grepl(look, anno$pred))

	if (length(x) == 0) next

	a = which(anno$SOterm[x] != "")
	b = which(anno$SOterm[x] == "")

	if (length(a) > 0) anno$SOterm[x[a]] = sprintf("%s|%s", anno$SOterm[x[a]], look)
	if (length(b) > 0) anno$SOterm[x[b]] = look

	if (length(a) > 0) anno$impact[x[a]] = sprintf("%s|%s", anno$impact[x[a]], imp)
	if (length(b) > 0) anno$impact[x[b]] = imp
}

del = which(anno$SOterm == "")
if (length(del) > 0) anno = anno[-del, ]


x = grepl("|", anno$impact)
if (length(x) > 0) {
	z = strsplit(anno$impact[x], "|", T)
	z = sapply(z, function(x) {
		x = unique(x)
		if (length(x) == 1)
			return(x[1])
		paste(x, collapse = "|")
	})
	anno$impact[x] = z
}


anno$delete = NA
anno$damage = NA
anno$benign = NA

x = which(grepl("deleterious", anno$pred))
if (length(x) > 0) anno$delete[x] = as.numeric(sub(".*deleterious\\(([\\.0-9]+)\\).*", "\\1", anno$pred[x]))

x = which(grepl("probably_damaging", anno$pred))
if (length(x) > 0) anno$damage[x] = as.numeric(sub(".*probably_damaging\\(([\\.0-9]+)\\).*", "\\1", anno$pred[x]))

x = which(grepl("benign", anno$pred))
if (length(x) > 0) anno$benign[x] = as.numeric(sub(".*benign\\(([\\.0-9]+)\\).*", "\\1", anno$pred[x]))


anno$EAS_AF = NA
anno$AMR_AF = NA
anno$AFR_AF = NA
anno$EUR_AF = NA
anno$SAS_AF = NA

x = which(grepl("EAS_AF", anno$pred))
if (length(x) > 0) anno$EAS_AF[x] = as.numeric(sub(".*EAS_AF=([\\.0-9]+).*", "\\1", anno$pred[x]))

x = which(grepl("AMR_AF", anno$pred))
if (length(x) > 0) anno$AMR_AF[x] = as.numeric(sub(".*AMR_AF=([\\.0-9]+).*", "\\1", anno$pred[x]))

x = which(grepl("AFR_AF", anno$pred))
if (length(x) > 0) anno$AFR_AF[x] = as.numeric(sub(".*AFR_AF=([\\.0-9]+).*", "\\1", anno$pred[x]))

x = which(grepl("EUR_AF", anno$pred))
if (length(x) > 0) anno$EUR_AF[x] = as.numeric(sub(".*EUR_AF=([\\.0-9]+).*", "\\1", anno$pred[x]))

x = which(grepl("SAS_AF", anno$pred))
if (length(x) > 0) anno$SAS_AF[x] = as.numeric(sub(".*SAS_AF=([\\.0-9]+).*", "\\1", anno$pred[x]))


save(anno, file = sprintf("raw.anno.chr%s.RData", CHR))

}



######################

library(data.table)


pack.data = function(d, chr, dir, grp, size = 50, nmax = 100) {
	
	for (tag in names(d)) {
		
		if (grepl("\\|", tag)) next
		
		out = d[[tag]]
		
		if (nrow(out) > 1) {
			
			out = out[sample(1:nrow(out)), ]
			
			if (nrow(out) < size * 1.5) {
				
				cat(sprintf(" %d  %s  %s  %d\n", chr, grp, tag, nrow(out)))
				cat(out$pos, file = sprintf("%s/%s.chr%d.%s.0.txt", dir, grp, chr, tag), sep = "\n")
				
			} else {
				
				spt = ceiling(nrow(out) / size)
				if (spt > nmax) {
					spt = nmax
				}
				out = split(out, cut(1:nrow(out), spt))
				
				for (sub in 1:length(out)) {
					
					if (nrow(out[[sub]]) > 100) {
						out[[sub]] = out[[sub]][1:100, ]
					}
					
					cat(sprintf(" %d  %s  %s  %d\n", chr, grp, tag, nrow(out[[sub]])))
					cat(out[[sub]]$pos, file = sprintf("%s/%s.chr%d.%s.%d.txt", dir, grp, chr, tag, sub-1), sep = "\n")
					
				}
			}
		}
		
	}
	
	save(d, file = sprintf("%s/_%s.chr%s.RData", dir, grp, chr))
}


for (chr in 1:22) {
	print(chr)

	chr.dir = sprintf("./chr%d", chr)
	dir.create(chr.dir, showWarnings = F)

	M = fread(sprintf("../data/1000G.chr%d.marker.txt", chr), header = T, stringsAsFactors = F)

	load(sprintf("raw.anno.chr%d.RData", chr))
	S = anno

	x = which(S$pos %in% M$Position)
	y = which(S$rsid %in% M$Label)
	z = intersect(x, y)

	if (length(z) > 0) {

		M = M[z, ]

		x = which(M$GenotypeCount1 > 1 & M$GenotypeCount1 < 51 & M$GenotypeCount2 == 0)

		if (length(x) > 0) {

			M = M[x, ]

			p = intersect(M$Position, S$pos)

			if (length(p) > 0) {

				x = which(S$pos %in% p)
				S = S[x, ]
				
				trm = split(S, S$SOterm)
				imp = split(S, S$impact)
				del = split(S, 100 * round(S$delete, 2))
				dam = split(S, 100 * round(S$damage, 2))
				ben = split(S, 100 * round(S$benign, 1))
				
				
				pack.data(trm, chr, chr.dir, "TRM", 10, 100)
				pack.data(imp, chr, chr.dir, "IMP", 10, 100)
				pack.data(del, chr, chr.dir, "DEL", 10)
				pack.data(dam, chr, chr.dir, "DAM", 10)
				pack.data(ben, chr, chr.dir, "BEN", 10)

			}
		}
	}
}











