


library(data.table)




bed = "20140520.strict_mask.autosomes.bed"

bed = fread(bed, header = F, stringsAsFactors = F)
bed$V4 = NULL

names(bed) = c("chr", "beg", "end")

bed = split(bed, bed$chr)





files = dir(pattern = "^ALL\\.chr([0-9]+)\\..+\\.vcf\\.gz$")

for (file in files) {
	
	chr = (sub("^ALL\\.chr([0-9]+)\\..+\\.vcf\\.gz$", "\\1", file))
	
	cat(chr, "\n")
	
	
	d = read.table(file, header = F, stringsAsFactors = F)
	d = as.data.table(d)
	
	names(d) = c("chr", "pos", "info", "a0", "a1", "qual", "fltr", "str")
	
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth at Site">
	##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Coverage, Less is worse">
	##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of Alleles in Samples with Coverage, Less is worse">
	##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Counts in Samples with Coverage">
	##INFO=<ID=AF,Number=.,Type=Float,Description="Alternate Allele Frequencies">
	##INFO=<ID=SVM,Number=1,Type=Float,Description="Frequency-adjusted SVM score">
	
	d$dp = as.numeric(sub("^DP=([0-9]+).+", "\\1", d$str))
	d$ns = as.numeric(sub(".+NS=([0-9]+).+", "\\1", d$str))
	d$an = as.numeric(sub(".+AN=([0-9]+).+", "\\1", d$str))
	d$ac = as.numeric(sub(".+AC=([0-9]+).+", "\\1", d$str))
	d$af = as.numeric(sub(".+AF=([0-9]+).+", "\\1", d$str))
	d$svm = as.numeric(sub(".+SVM=(.+)$", "\\1", d$str))
	
	
	m = fread(sprintf("../data/1000G.chr%s.marker.txt", chr), header = T, stringsAsFactors = F)
	
	if (any(duplicated(m$Position))) {
		del = which(duplicated(m$Position))
		m = m[-del, ]
	}
	
	if (any(duplicated(d$pos))) {
		del = which(duplicated(d$pos))
		d = d[-del, ]
	}
	
	
	key = intersect(m$Position, d$pos)
	if (length(key) == 0) next
	
	i = which(d$pos %in% key)
	if (length(i) == 0) next
	d = d[i,]
	
	i = which(m$Position %in% key)
	if (length(i) == 0) next
	m = m[i,]
	
	m = m[order(m$Position),]
	d = d[order(d$pos),]
	
	if (!identical(d$pos, m$Position)) stop("Mismatch!")
	
	
	p = cbind(d, m)
	
	
	x = which(p$AlleleCount1 > 1 & p$AlleleCount1 < 5000-500)
	if (length(x) == 0) stop("NONE!")
	p = p[x,]
	
	
	# x = which(p$ac > 1 & p$ac <= 5000)
	# if (length(x) == 0) stop("NONE!")
	# p = p[x,]
	
	
	x = which(p$dp > mean(d$dp))
	if (length(x) == 0) stop("NONE!")
	p = p[x,]
	
	
	x = which(p$ns > mean(d$ns))
	if (length(x) == 0) stop("NONE!")
	p = p[x,]
	
	
	# x = which(p$svm > mean(d$svm))
	# if (length(x) == 0) stop("NONE!")
	# p = p[x,]
	
	
	# x = which(p$ac >= p$AlleleCount1)
	# if (length(x) == 0) stop("NONE!")
	# p = p[x,]
	
	
	
	
	p$mask = sapply(p$pos, function(x, z) {
		length(which(z$beg < x & z$end > x))
	}, bed[[sprintf("chr%s", chr)]])
	
	if (length(which(p$mask == 1)) == 0) stop("NONE!")
	
	del = which(p$mask != 1)
	if (length(del) > 1) {
		p = p[-del,]
	}
	
	cat(sprintf("Sites: %d total, %d selected\n", nrow(d), nrow(p)))
	
	
	cat(p$pos, sep = "\n", file = sprintf("pos.chr%s.filtered.txt", chr))
	save(p, file = sprintf("pos.chr%s.filtered.RData", chr))
	
	
	# tmp = split(p, p$ac)
	# tmp = lapply(tmp, function(x) {
	# 	if (nrow(x) < 1000) return(x)
	# 	x[sample(1:nrow(x), 1000),]
	# })
	# tmp = rbindlist(tmp)
	# 
	# if (nrow(tmp) > 5000) {
	# 	tmp = tmp[sample(1:nrow(tmp), 5000),]
	# }
	# 
	# cat(sort(tmp$pos), sep = "\n", file = sprintf("expos.chr%s.n5000.txt", chr))
	
	
	cat("\n")
}

stop()



# packs

npp = 1000

files = dir(pattern = "pos.chr([0-9]+).filtered.txt")

for (file in files) {
	print(file)
	
	chr = as.numeric(sub("pos.chr([0-9]+).filtered.txt", "\\1", file))
	
	dir.create(sprintf("chr%d", chr), showWarnings = F)
	
	pos = as.numeric(readLines(file))
	key = cut(1:length(pos), breaks = ceiling(length(pos) / npp), include.lowest = T)
	tmp = split(pos, key)
	
	for (i in 1:length(tmp)) {
		cat(sort(tmp[[i]]), sep = "\n", file = sprintf("chr%d/pack%03d.chr%d.filtered.txt", chr, i, chr))
	}
}



stop()


for (file in dir(pattern = "^expos\\..+\\.RData$")) {
	tag = sub("^(.+)\\.RData$", "\\1", file)
	print(tag)
	
	load(file)
	
	x = which(p$ac >= 15)
	pos = p$pos[x]
	
	if (length(pos) > 10) {
		pos = sample(pos, size = 10)
	}
	
	pos = sort(pos)
	
	cat(pos, sep = "\n", file = sprintf("_n10.%s.txt", tag))
}


stop()









library(ggplot2)
library(ggthemes)

Ne = 2 * 10000


for (chr in 15:15) {
	cat("Chr:", chr, "")
	
	marker = fread(sprintf("../data/1000G.chr%d.marker.txt", chr), header = T)
	
	for (near in c("NN", "RD")) {
		
		site = fread(sprintf("__run_1kg_chr%d_%s_1000.cle.txt", chr, near), header = T)
		
		pair = fread(sprintf("__run_1kg_chr%d_%s_1000.ccf.txt", chr, near), header = T)
		
		pair = split(pair, list(pair$MarkerID))
		
		for (tag in names(pair)) {
			cat(".")
			
			x = pair[[ tag ]]
			
			m = marker[which(marker$MarkerID == x$MarkerID[1]),]
			title = sprintf("Chromosome = %d     Position = %s     ID = %s     Alleles = %s     Allele count = %d", m$Chromosome, m$Position, m$Label, m$Alleles, m$AlleleCount1)
			
			
			y = which(site$MarkerID == x$MarkerID[1])
			y = site[y,]
			
			y$fltr = "Raw"
			y$fltr[which(y$Adjusted == 1)] = "Adjusted"
			y$fltr = factor(y$fltr, levels = c("Raw", "Adjusted"), labels = c("Raw", "Adjusted"), ordered = T)
			
			x = rbindlist(lapply(split(x, x$Clock), function(x) {
				x = x[order(abs(x$Shared-1), x$q50),]
				
				a = which(x$Shared == 1)
				b = which(x$Shared == 0)
				
				x$i = NA
				x$i[a] = ((1:length(a)) / (length(a)+1)) - 0.06
				x$i[b] = 1 + ((1:length(b)) / (length(b)+1)) + 0.06
				
				x
			}))
			
			x$cord = NA
			x$cord[ which(x$Shared == 1) ] = "Concordant pairs"
			x$cord[ which(x$Shared == 0) ] = "Discordant pairs"
			
			x$Clock[ which(x$Clock == "M") ] = "(a) Mutation clock"
			x$Clock[ which(x$Clock == "R") ] = "(b) Recombination clock"
			x$Clock[ which(x$Clock == "C") ] = "(c) Combined clock"
			
			y$Clock[ which(y$Clock == "M") ] = "(a) Mutation clock"
			y$Clock[ which(y$Clock == "R") ] = "(b) Recombination clock"
			y$Clock[ which(y$Clock == "C") ] = "(c) Combined clock"
			
			
			# q = split(x, x$Clock)
			# q = lapply(q, function(x) {
			# 	s = which(x$Shared == 1)
			# 	u = which(x$Shared == 0)
			# 
			# 	s = x[s,]
			# 	u = x[u,]
			# 
			# 	# del = which(u$SegmentRHS - u$SegmentLHS < 5)
			# 	# if (length(del) > 0) {
			# 	# 	u = u[-del,]
			# 	# }
			# 
			# 	s$con0 = sprintf("%d.%d", s$SampleID0, s$Chr0)
			# 	s$con1 = sprintf("%d.%d", s$SampleID1, s$Chr1)
			# 	u$dis = sprintf("%d.%d", u$SampleID0, u$Chr0)
			# 	
			# 	uni = unique(c(s$con0, s$con1, u$dis))
			# 	out = NULL
			# 	
			# 	for (i in 1:length(uni)) {
			# 		z = seq(0, 2, length.out = length(uni))[i]
			# 		a = which(s$con0 == uni[i])
			# 		b = which(s$con1 == uni[i])
			# 		c = which(u$dis == uni[i])
			# 		if (length(a) > 0) out = rbind(out, data.table(Clock = x$Clock[1], x0 = z, x1 = s$i[a], y0 = 1e5, y1 = s$q50[a]))
			# 		if (length(b) > 0) out = rbind(out, data.table(Clock = x$Clock[1], x0 = z, x1 = s$i[b], y0 = 1e5, y1 = s$q50[b]))
			# 		if (length(c) > 0) out = rbind(out, data.table(Clock = x$Clock[1], x0 = z, x1 = u$i[c], y0 = 1e5, y1 = u$q50[c]))
			# 	}
			# 
			# 	# for (i in 1:nrow(u)) {
			# 	# 	z = unique(which(s$con0 == u$dis[i]), which(s$con1 == u$dis[i]))
			# 	# 	if (length(z) == 0) next
			# 	# 	out = rbind(out, data.table(Clock = x$Clock[1], x0 = s$i[z], x1 = u$i[i], y0 = s$q50[z], y1 = u$q50[i]))
			# 	# }
			# 	out
			# })
			# q = rbindlist(q)
			
			
			q = split(x, list(x$Shared, x$Clock))
			q = lapply(q, function(x) {
				i = which(x$Pass == 0)
				if (length(i) == 0) return(data.table(Clock = x$Clock[1], cord = x$cord[1], x = -10))
				if (x$Shared[1] == 0) {
					i = max(i)
					if (i == nrow(x)) return(data.table(Clock = x$Clock[1], cord = x$cord[1], x = -10))
					return(data.table(Clock = x$Clock[1], cord = x$cord[1], x = x$i[i+1] + abs(x$i[i+1] - x$i[i])))
				}
				if (x$Shared[1] == 1) {
					i = min(i)
					if (i == 1) return(data.table(Clock = x$Clock[1], cord = x$cord[1], x = -10))
					return(data.table(Clock = x$Clock[1], cord = x$cord[1], x = x$i[i-1] + abs(x$i[i-1] - x$i[i])))
				}
			})
			q = rbindlist(q)
			
			
			tmp = split(x, x$Clock)
			tmp = lapply(tmp, function(z) {
				z$excl = F
				a = z[which(z$Shared == 1),]
				b = z[which(z$Shared == 0),]
				for (i in 1:nrow(b)) {
					j = which(a$SampleID0 == b$SampleID0[i] | a$SampleID1 == b$SampleID0[i])
					if (length(j) == 0) next
					k = which(a$q50[j] > b$q50[i])
					if (length(k) > 0) {
						b$excl[i] = T
						a$excl[j[k]] = T
					}
				}
				rbind(a, b)
			})
			tmp = rbindlist(tmp)
			del = which(tmp$excl)
			if (length(del) > 0) {
				tmp = tmp[-del,]
			}
			
			
			gg = ggplot(x) +
				facet_grid(.~Clock) +
				geom_rect(data = data.table(xmin = 0.95, xmax = 1.05, ymin = 1e-8, ymax = 1e8), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", alpha = 1/2) +
				geom_ribbon(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, fill = cord), alpha = 1/2, size = 1/3, show.legend = F) +
				geom_line(aes(x=i, y=q50 * Ne, color = cord), size = 1/3, alpha = 1, show.legend = F) +
				#geom_linerange(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, color = concord), alpha = 1/2, size = 1/3, show.legend = F) +
				#geom_point(aes(x=i, y=q50 * Ne, color = concord), shape = 16, alpha = 1, size = 1/2, show.legend = F) +
				#geom_vline(data = q, aes(xintercept = x), size = 1/3, alpha = 1/2) +
				geom_rect(data = y, aes(xmin = 0.955, xmax = 1.045, ymin=PostCI025 * Ne, ymax=PostCI975 * Ne, fill = fltr), alpha = 1/2) +
				geom_segment(data = y, aes(x = 0.956, xend = 1.044, y = PostMode * Ne, yend = PostMode * Ne, linetype = fltr), alpha = 1/2, size = 1/2) +
				#geom_segment(data = q, aes(x = x0, xend = x1, y = y0, yend = y1 * Ne), alpha = 1/10, size = 1/4) +
				scale_x_continuous(breaks = c(0.45, 1, 1.55), labels = c("Concordant\npairs", "Age\nestimate", "Discordant\npairs")) +
				scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ','), minor_breaks = c(seq(1,10,by=1), seq(10,100,by=10), seq(100,1000,by=100), seq(1000,10000,by=1000), seq(10000,100000,by=10000))) +
				scale_colour_manual(values = c("Concordant pairs" = "forestgreen", "Discordant pairs" = "chocolate", "Raw" = "goldenrod", "Adjusted" = "dodgerblue")) +
				scale_fill_manual(values = c("Concordant pairs" = "forestgreen", "Discordant pairs" = "chocolate", "Raw" = "goldenrod", "Adjusted" = "dodgerblue")) +
				scale_linetype_manual(values = c("Raw" = "11", "Adjusted" = "solid")) +
				coord_cartesian(ylim = c(0.8, 1.2e5), xlim = c(-0.075, 2.075), expand = F) +
				theme_few() +
				theme(aspect.ratio = 16/9,
							panel.background = element_rect(fill = "grey90"),
							panel.border = element_rect(fill = NA, colour = "black", size = 2/3),
							strip.text.x = element_text(face = "bold", hjust = 0, size = 12),
							axis.title.x = element_blank(),
							axis.ticks.x = element_blank(),
							axis.text.x = element_text(size = 8),
							axis.title.y = element_text(margin = margin(1,-1,1,1, "mm")),
							plot.title = element_text(size = 14),
							plot.margin = margin(1, 1, 1, 1, "cm"),
							axis.line = element_blank(),
							legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
							legend.key.size = unit(4.25, "mm"), 
							legend.text = element_text(size = 8),
							legend.justification = c(1, 0), legend.position = c(0.995, 1-0.991), 
							legend.title = element_blank(),
							legend.margin = margin(-0.25,2,1,1, "mm"),
							panel.grid.major.y = element_line(colour = "grey99", size = 1/2),
							panel.grid.minor.y = element_line(colour = "grey99", size = 1/3),
							panel.grid.major.x = element_blank(),
							panel.grid.minor.x = element_blank()) +
				ylab("Time (generations)") +
				ggtitle(title)
			gg
			
			ggsave(gg, filename = sprintf("_plot_run_ccf.1kg_chr%d.%d.%s.pdf", chr, x$MarkerID[1], near), width = 297, height = 210, units = "mm")
			
			
		}
		
	}
	
	cat("\n")
}










files = dir("^pack.+", path = "./chr20", full.names = T)
pos = c()
for (file in files) {
	pos = c(pos, as.numeric(readLines(file)))
}


i = match(pos, marker$Position)


m = marker[i,]

k = which(m$AlleleCount1 <= 50)

m = m[k,]


m = m[order(sample(1:nrow(m))),]

z = split(m$Position, cut(1:nrow(m), 12))
lapply(z, function(x) {
	x = sort(x)
	cat(x, file = sprintf("ibd_pack.%s.txt", paste(sample(c(0:9, letters, LETTERS), 16), collapse = "")), sep = "\n")
})







