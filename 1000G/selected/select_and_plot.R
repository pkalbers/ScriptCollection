
rs = c(
"rs1874165",
"rs2973631",
"rs3805547",
"rs6875787",
"rs16874441",
"rs55862350",
"rs56001636",
"rs56256550",
"rs58945509",
"rs58979818",
"rs61051796",
"rs1994929",
"rs2099081",
"rs2115262",
"rs2914281")

rs = unique(c(
"rs7714612",
"rs10042274",
"rs113958470",
"rs141880427",
"rs181064967",
"rs182972078",
"rs183034145",
"rs184600328",
"rs185548599",
"rs186121405",
"rs188517683",
"rs200723191",
"rs186504490",
"rs146505774"))


pos = c()

for (r in rs) {
	i = which(marker$Label == r)
	if (length(i) != 0) {
		cat(r, "\n")
		print(marker[i,])
		cat("\n")
		pos = c(pos, marker$Position[i])
	}
}

cat(pos, file = "pos.PRDM9.more.txt", sep = "\n")











library(ggplot2)
library(ggthemes)

Ne = 2 * 10000


chr = 5


marker = fread(sprintf("../data/1000G.chr%d.marker.txt", chr), header = T)



for (iter in 2:3) {
for (near in c("NN", "RD")) {


ccf = fread(sprintf("PRDM9_more_1kg_chr%d.%d.%s.ccf.txt", chr, iter, near), header = T)
cle = fread(sprintf("PRDM9_more_1kg_chr%d.%d.%s.cle.txt", chr, iter, near), header = T)


ccf = split(ccf, list(ccf$MarkerID))


for (tag in names(ccf)) {
	x = ccf[[ tag ]]
	
	
	m = marker[which(marker$MarkerID == x$MarkerID[1]),]
	title = sprintf("Chromosome = %d     Position = %s     ID = %s     Alleles = %s     Allele count = %d", m$Chromosome, m$Position, m$Label, m$Alleles, m$AlleleCount1)
	
	
	y = which(cle$MarkerID == x$MarkerID[1])
	y = cle[y,]
	
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
	
	
	gg = ggplot(x) +
		facet_grid(.~Clock) +
		geom_rect(data = data.table(xmin = 0.95, xmax = 1.05, ymin = 1e-8, ymax = 1e8), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", alpha = 1/2) +
		geom_ribbon(data = x[which(x$Shared == 1),], aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne), fill = "forestgreen", alpha = 1/2, size = 1/3, show.legend = F) +
		geom_line(data = x[which(x$Shared == 1),], aes(x=i, y=q50 * Ne), color = "forestgreen", size = 1/3, alpha = 1, show.legend = F) +
		geom_ribbon(data = x[which(x$Shared == 0),], aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne), fill = "chocolate", alpha = 1/2, size = 1/3, show.legend = F) +
		geom_line(data = x[which(x$Shared == 0),], aes(x=i, y=q50 * Ne), color = "chocolate", size = 1/3, alpha = 1, show.legend = F) +
		#geom_linerange(aes(x=i, ymin=q25 * Ne, ymax=q75 * Ne, color = cord), alpha = 1/2, size = 1/3, show.legend = F) +
		#geom_point(aes(x=i, y=q50 * Ne, color = cord), shape = 16, alpha = 1, size = 1/2, show.legend = F) +
		geom_rect(data = y, aes(xmin = 0.955, xmax = 1.045, ymin=PostCI025 * Ne, ymax=PostCI975 * Ne, fill = fltr), alpha = 2/3) +
		geom_segment(data = y, aes(x = 0.956, xend = 1.044, y = PostMode * Ne, yend = PostMode * Ne, linetype = fltr), alpha = 1, size = 1/2) +
		scale_x_continuous(breaks = c(0.45, 1, 1.55), labels = c("Concordant\npairs", "Age\nestimate", "Discordant\npairs")) +
		scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), scientific = F, big.mark = ','), minor_breaks = c(seq(1,10,by=1), seq(10,100,by=10), seq(100,1000,by=100), seq(1000,10000,by=1000), seq(10000,100000,by=10000))) +
		scale_fill_manual(values = c("Raw" = "goldenrod", "Adjusted" = "dodgerblue")) +
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
	
	
	ggsave(filename = sprintf("__plot_PRDM9.chr%d.%d.%s.%d.pdf", chr, x$MarkerID[1], near, iter), width = 297, height = 210, units = "mm")
	
}

}}




