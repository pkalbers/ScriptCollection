

library(data.table)
library(ggplot2)
library(ggthemes)


Ne = 2 * 7300 # 10000

marker = fread("../truH.marker.txt", header = T)
times = fread("~/Research/africa/dev/OutOfAfricaHapMap20.times.txt", header = T)



rmsle = function(a, b) {
	n = length(a)
	a = log10(a)
	b = log10(b)
	sqrt(sum((a - b)^2) / n)
}

mle = function(a, b) {
	a = log10(a)
	b = log10(b)
	10^(mean(abs(a - b)))
}

se <- function(x) sqrt(var(x)/length(x))

expected.age = function(k, n) {
	x = k / n
	if (x == 0) return(0)
	if (x == 1) return(2)
	((-2 * x) / (1 - x)) * log(x)
}



d = NULL

for (file in dir(pattern = "\\.cle\\.txt$", path = ".", full.names = T)) {
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	tmp = cbind(tmp, times[tmp$MarkerID + 1, ])
	
	tag = sub("^ibd_con_dis\\.(.+)\\.cle.+$", "\\1", basename(file))

	
	tmp$data = ""
	tmp$phased = ""
	#tmp$expos = ""
	#tmp$tmrca = as.numeric(sub("^.+_([0-9]+)_.+$", "\\1", tag))
	tmp$near = ""
	#tmp$fltr = ""
	
	if (grepl("^tru_.+$", tag)) tmp$data = "(a) Before error"
	if (grepl("^err_.+$", tag)) tmp$data = "(b) After error"
	
	if (grepl("^.+_H.+$", tag)) tmp$phased = " Simulated haplotypes "
	if (grepl("^.+_P.+$", tag)) tmp$phased = "Phased haplotypes"
	
	# if (grepl("^.+_id_.+$", tag)) tmp$expos = "(a) Identical"
	# if (grepl("^.+_fp_.+$", tag)) tmp$expos = "(b) False positives"
	# if (grepl("^.+_fn_.+$", tag)) tmp$expos = "(c) False negatives"
	# if (grepl("^.+_ff_.+$", tag)) tmp$expos = "(d) Both"
	
	if (grepl("_NN", tag)) tmp$near = "Nearest neighbour"
	if (grepl("_RD", tag)) tmp$near = "Random selection"
	
	#if (grepl("_none$", tag)) tmp$fltr = " Raw "
	#if (grepl("_prop$", tag)) tmp$fltr = "(b) Fixed proportion"
	#if (grepl("_auto$", tag)) tmp$fltr = "Adjusted"

	tmp$Clock[which(tmp$Clock == "M")] = "(a) Mutation clock"
	tmp$Clock[which(tmp$Clock == "R")] = "(b) Recombination clock"
	tmp$Clock[which(tmp$Clock == "C")] = "(c) Combined clock"
	
	tmp$tag = tag
	
	d = rbind(d, tmp)
}


p = d

p = p[(which(p$data == "(a) Before error")), ]
p = p[(which(p$data == "(b) After error")), ]

p = p[(which(p$Adjusted == 0)), ]
p = p[(which(p$Adjusted == 1)), ]

p = p[(which(p$near == "Nearest neighbour")), ]
p = p[(which(p$near == "Random selection")), ]




D = split(d, list(d$near))



for (tag in names(D)) {
print(tag)
break
p = D[[tag]]

#p = p[which(p$Fk <= 25),]
#p = p[which(p$N_Others >= 50),]
#p = p[which(p$N_Shared >= choose(p$Fk, 2)),]


p$exp = sapply(p$Fk, expected.age, 5000) * 10000 * 2
p$mid = exp(log(p$node.time + 1e-8) + ((log(p$parent.time + 1e-8) - log(p$node.time + 1e-8)) / 2))
p$x = p$mid # p$node.time


# z = which(p$PostMode * Ne >= p$node.time & p$PostMode * Ne <= p$parent.time)
# p$x[z] = p$PostMode[z] * Ne
# 
# z = which(p$PostMode * Ne < p$node.time)
# p$x[z] = p$node.time[z]
# 
# z = which(p$PostMode * Ne > p$parent.time)
# p$x[z] = p$parent.time[z]


mx = c(0.5, 8e4)

tm = rbind(data.table(x = 10^(0:6), y = 1e-6, xend = 10^(0:6), yend = 0.7),
					 data.table(y = 10^(0:6), x = 1e-6, yend = 10^(0:6), xend = 0.7))

tn = c((2:9), 10*(2:9), 100*(2:9), 1000*(2:9), 10000*(2:9))
tn = rbind(data.table(x = tn, xend = tn, y = 1e-6, yend = 0.615), 
					 data.table(y = tn, yend = tn, x = 1e-6, xend = 0.615))

tx = rbind(data.table(x = 1e-6, xend = 1e6, y = 1e-6, yend = 0.7),
					 data.table(x = 1e-6, xend = 0.7, y = 1e-6, yend = 1e6))


w = split(p, sprintf("%s %s", p$Clock, p$tag))
w = lapply(w, function(x) {
	x$rankc = cor(x$x, x$PostMode, method = "spearman")
	x$rmsle = rmsle(x$x + 1e-32, x$PostMode * Ne + 1e-32)
	x[1,]
})
w = rbindlist(w)

w$rankc.txt = sprintf("Spearman's rho = %.3f", w$rankc)
w$rmsle.txt = sprintf("RMSLE = %.3f", w$rmsle)



scat = ggplot(p) +
	facet_grid(phased~Clock) +
	geom_raster(aes(x, PostMode * Ne, fill = (..density.. / max(..density..)) * 100 ), position = "identity", stat = "bin2d", binwidth = c(0.075, 0.075)) +
	geom_tile(aes(x, PostMode * Ne, fill = (..density.. / max(..density..)) * 100), position = "identity", stat = "bin2d", binwidth = c(0.075, 0.075)) +
	geom_smooth(aes(x, node.time + 1e-8), method = "lm", size = 1/2, color = "grey30", fullrange = F, se = T) +
	#geom_smooth(aes(x, mid), method = "lm", size = 1/2, color = "grey20", fullrange = T) +
	geom_smooth(aes(x, parent.time + 1e-8), method = "lm", size = 1/2, color = "grey30", fullrange = F, se = T) +
	geom_text(data = w, aes(x = 50000, y = 2.4, label = rankc.txt), size = 3, hjust = 1) +
	geom_text(data = w, aes(x = 50000, y = 1.2, label = rmsle.txt), size = 3, hjust = 1) +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = F) +
	geom_smooth(aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", se = T, fullrange = F) +
	geom_rect(data = tx, aes(xmin=x, xmax=xend, ymin=y, ymax=yend), fill="grey80") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), size = 1/2, color="black") +
	geom_segment(data = tn, aes(x=x, xend=xend, y=y, yend=yend), size = 1/4, color="black") +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), na.value = "grey90", trans = "log10", limits = c(1, 100), breaks = c(1, 10, 100), labels = c("<1%", "10%", "100%")) +
	scale_color_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy"), na.value = "grey90", trans = "log10", limits = c(1, 100), breaks = c(1, 10, 100), labels = c("<1%", "10%", "100%")) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 9),
				#axis.ticks = element_blank(),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(size = 10),
				legend.justification=c(0,1), legend.position=c(1-0.98,0.99),
				legend.background = element_rect(fill = NA, colour = NA, size = NA),
				#legend.position = "bottom",
				plot.title = element_text(face = "bold", size = 12),
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 7),
				legend.title = element_blank()) +
	ylab("Inferred age (generations)") + xlab("True age (generations)") +
	ggtitle(p$data[1])
scat

ggsave(scat, filename = ("__plot.cle.tru_NN.pdf"), height = 7.5, width = 11)
ggsave(scat, filename = ("__plot.cle.tru_RD.pdf"), height = 7.5, width = 11)

ggsave(scat, filename = ("__plot.cle.err_NN.pdf"), height = 7.5, width = 11)
ggsave(scat, filename = ("__plot.cle.err_RD.pdf"), height = 7.5, width = 11)


ggsave(scat, filename = sprintf("_B_plot_sites.scat.%s.pdf", make.names(p$tag[1])), height = 5.5, width = 8)

}

p$intv = sprintf("           %s", (as.character(cut(p$Fk / 5000 * 100, breaks = (0:10)/10, labels = (1:10)/10, include.lowest = T))))
p$intv = sprintf("          %d", as.numeric(as.character(cut(p$node.time, breaks = 10^(0:5), labels = 10^(1:5), include.lowest = T))))
p$intv = cut(p$node.time, breaks = 10^seq(log10(1), log10(1e5), length.out = 26), include.lowest = T)

p$est = NA

z = which(p$PostMode * Ne >= p$node.time & p$PostMode * Ne <= p$parent.time &  p$PostMode * Ne >= p$x)
p$est[z] = 1 + 1e-8
z = which(p$PostMode * Ne >= p$node.time & p$PostMode * Ne <= p$parent.time &  p$PostMode * Ne < p$x)
p$est[z] = 1 - 1e-8

z = which(p$PostMode * Ne > p$parent.time)
p$est[z] = p$PostMode[z] * Ne / p$parent.time[z]

length(z) / nrow(p) * 100

z = which(p$PostMode * Ne < p$node.time)
p$est[z] = p$PostMode[z] * Ne / p$node.time[z]

length(z) / nrow(p) * 100

vl = data.table(x = (2:26) - 0.5)
hl = NULL
for (y in c(1/(20:2), 1, (2:20))) {
	hl = rbind(hl, data.table(y = y, x0 = (1:11) - 0.6, x1 = (1:11) - 0.4))
}

del = which(is.na(p$intv))
if (length(del) > 0) {
	p = p[-del,]
}

ggplot(p) +
	facet_grid(data~Clock) +
	geom_raster(aes(intv, est, fill = ..density.. * 100), position = "identity", stat = "bin2d", binwidth = c(1/26, 1/20)) +
	#geom_tile(aes(intv, est, fill = ..density.. * 100), position = "identity", stat = "bin2d", binwidth = c(10/11, 1/25)) +
	#geom_violin(aes(intv, est), scale = "area", draw_quantiles = c(0.25, 0.5, 0.75), trim = F, color = "grey20", size = 1/3) +
	geom_vline(data = vl, aes(xintercept = x), size = 1/3, color = "grey50", alpha = 1/4) +
	geom_hline(data = hl, aes(yintercept = y), size = 1/3, color = "grey50", alpha = 1/4) +
	#geom_segment(data = hl, aes(x = x0, xend = x1, y = y, yend = y), size = 1/3, color = "grey50", alpha = 1/2) +
	coord_cartesian(ylim = c(0.05, 20), xlim = c(0.5,22.5), expand = F) +
	#scale_x_log10() +
	#scale_y_log10(breaks = 10^(-2:2), labels = trimws(format(10^(-2:2), big.mark = ',', scientific = F, drop0trailing = T))) +
	scale_y_log10(breaks = c(1/10, 1/(5:2), 1, (2:5), 10), labels = format(round(c(1/10, 1/(5:2), 1, (2:5), 10), 3), trim = T, scientific = F, drop0trailing = T)) +
	#scale_fill_gradientn(colours = c("grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy", "black"), na.value = "grey90",
	#										 trans = "log10", limits = c(0.1, 100), breaks = c(0.1,1,10,100), labels = c("0.1%", "1%", "10%", "100%")) +
	scale_fill_gradientn(colours = c("grey90", "khaki1", "khaki2", "orange", "orangered", "purple2", "purple4", "black"), na.value = "grey90", 
											 trans = "log10", limits = c(0.1, 100), breaks = c(0.1,1,10,100), labels = c("0.1%", "1%", "10%", "100%")) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 8),
				#axis.ticks = element_blank(),
				panel.background = element_rect(fill = "grey90"),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				axis.ticks.x = element_blank(),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,1), legend.position=c(1-0.99333,0.99), legend.direction = "horizontal",
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
				#legend.position = "bottom",
				legend.key.height = unit(0.333, "cm"),
				legend.key.width = unit(0.75, "cm"),
				legend.margin = margin(1,3,1,2, "mm"),
				legend.text = element_text(size = 7.5),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("Allele frequency bin (%)")

ggsave(file = "__plot.pdf", height = 5.5, width = 8)



p$intv = p$parent.time - p$node.time

ggplot(p) +
	facet_grid(data~Clock) +
	geom_segment(aes(x = node.time, xend = parent.time, y = PostMode * Ne, yend = PostMode * Ne, group = MarkerID), alpha = 1/50) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) 


by(p, list(p$data, p$Clock), function(x) { cor(x$intv, x$PostMode, method = 's') })


p$intv = as.numeric(cut(p$Fk / 5000 * 100, breaks = (0:10)/10, labels = (1:10)/10, include.lowest = T))


a = lapply(split(p, list(p$intv, p$fltr, p$near, p$Clock)), function(x) {
	if (nrow(x) == 0) return(NULL)
	rbind(data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = "  tc  ", y = cor(x$PostMode * Ne, x$node.time, method = "s")),
				data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = " tm ", y = cor(x$PostMode * Ne, x$mid, method = "s")),
				data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = "td", y = cor(x$PostMode * Ne, x$parent.time, method = "s")))
})
a = rbindlist(a)

aa = ggplot(a) +
	facet_grid(type~clock, drop = T) +
	geom_line(aes(intv, y, color = factor(tag1), linetype = factor(tag2)), alpha = 0.8) +
	geom_point(aes(intv, y, color = factor(tag1), shape = factor(tag2)), alpha = 0.8) +
	coord_cartesian(ylim=c(-0.01, 1.01), xlim = c(0.5, 10.5), expand = F) +
	scale_x_continuous(breaks = (1:10), labels = (1:10)/10) +
	scale_y_continuous(breaks = (1:9)/10) +
	theme_bw() +
	theme(#aspect.ratio = 2/3,
				axis.text = element_text(size = 8),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", size = 9),
				strip.text.y = element_text(face = "bold", angle = -90, size = 9),
				#legend.justification=c(1,0), legend.position=c(0.99,1-0.99),
				#legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
				legend.position = "bottom",
				legend.text = element_text(size = 9),
				legend.title = element_blank()) +
	xlab("Allele frequency bin (%)") + ylab("Spearman's rho")
aa

ggsave(aa, filename = sprintf("_plot.stat_rankcor.%s.pdf", make.names(tag)), height = 7.5, width = 7.5)



b = lapply(split(p, list(p$intv, p$fltr, p$near, p$Clock)), function(x) {
	if (nrow(x) == 0) return(NULL)
	rbind(data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = "  tc  ", y = rmsle(x$PostMode * Ne, x$node.time)),
				data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = " tm ", y = rmsle(x$PostMode * Ne, x$mid)),
				data.table(intv = x$intv[1], tag1 = x$fltr[1], tag2 = x$near[1], clock = x$Clock[1], type = "td", y = rmsle(x$PostMode * Ne, x$parent.time)))
})
b = rbindlist(b)

bb = ggplot(b) +
	facet_grid(type~clock, drop = T) +
	geom_line(aes(intv, y, color = factor(tag1), linetype = factor(tag2)), alpha = 0.8) +
	geom_point(aes(intv, y, color = factor(tag1), shape = factor(tag2)), alpha = 0.8) +
	coord_cartesian(ylim=c(-0.01, 1.95), xlim = c(0.5, 10.5), expand = F) +
	scale_x_continuous(breaks = (1:10), labels = (1:10)/10) +
	scale_y_continuous(breaks = (1:20)/5) +
	theme_bw() +
	theme(#aspect.ratio = 2/3,
				axis.text = element_text(size = 8),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", size = 9),
				strip.text.y = element_text(face = "bold", angle = -90, size = 9),
				#legend.justification=c(1,1), legend.position=c(0.99,0.99),
				#legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
				legend.position = "bottom",
				legend.text = element_text(size = 9),
				legend.title = element_blank()) +
	xlab("Allele frequency bin (%)") + ylab("RMSLE")
bb

ggsave(bb, filename = sprintf("_plot.stat_rmsle.%s.pdf", make.names(tag)), height = 7.5, width = 7.5)





p$intv = (cut(p$Fk / 5000 * 100, breaks = (0:10)/10, include.lowest = T))

q = data.table(intv = p$intv, clock = p$Clock, tag1 = p$fltr, tag2 = p$near, est = p$PostMode * Ne, tru = p$mid, t0 = p$node.time, t1 = p$parent.time)

q$map = log(q$est / q$t0) / log(q$t1 / q$t0)

tmp = q
#tmp$fkr = cut(tmp$fk, breaks = c(2, 5, 10, 20, 30, 40, 50), include.lowest = T)
tmp = split(tmp, list(tmp$clock, tmp$tag1, tmp$tag2, tmp$intv))
tmp = lapply(tmp, function(z) {
	if (nrow(z) == 0) return(NULL)
	x = seq(-2, 3, by=0.05)
	y = ecdf(z$map)(x)
	data.table(intv = z$intv[1], clock = z$clock[1], tag1 = z$tag1[1], tag2 = z$tag2[1], x=x, y=y)
})
q = rbindlist(tmp)


hist = ggplot(q) + 
	facet_grid(tag2+tag1~clock) +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	geom_vline(xintercept = c(0,1), size = 2/4, linetype = "21") +
	geom_line(aes(x, y, color = factor(intv)), size = 1, alpha = 3/4) +
	geom_hline(yintercept = c(0,1), color = "grey60") +
	#geom_line(data = r, aes(x, y), color = "black", size = 1, alpha = 1/3) +
	scale_x_continuous(breaks = seq(-10, 10, by = 0.5), labels = as.character(seq(-10, 10, by = 0.5))) +
	scale_y_continuous(breaks = (1:9)/10) +
	#scale_y_continuous(breaks = seq(0, 10000, by = 250), labels = format(seq(0, 10000, by = 250), big.mark = ",")) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey60", "maroon1"))(3)[2],
																										 colorRampPalette(c("grey70", "navy"))(3)[2],
																										 colorRampPalette(c("grey50", "red"))(3)[2],
																										 colorRampPalette(c("grey95", "darkorange3"))(3)[2])))(length(unique(q$intv)))) +
	coord_cartesian(expand = F, xlim = c(-1.75, 2.75), ylim = c(-0.005, 1.005)) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title = element_text(size = 6.5),
				legend.text = element_text(size = 6),
				legend.key.height = unit(2.5, "mm"),
				legend.key.width  = unit(2.5, "mm"),
				legend.box.background = element_rect(fill = "white", size = 1/3, colour = "grey20"),
				legend.position = c(0.005, 0.997), legend.justification = c(c(0,1)), #legend.margin = margin(0,1,1,1, "mm"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.y = element_blank()) + 
	xlab("Relative age estimate (log-scale)") + ylab("CDF") +
	guides(color = guide_legend(title = "AF bin (%)", override.aes = list(size = 3, alpha = 1)))
hist

ggsave(hist, filename = sprintf("_plot.hist.%s.pdf", make.names(tag)), height = 12.5, width = 7.5)




# p$intv = (cut(p$Fk / 5000 * 100, breaks = (0:10)/10, include.lowest = T))
# p$conv = (p$PostCI975 * Ne) - (p$PostCI025 * Ne)
# 
# #del = which(p$conv < 1)
# #if (length(del) > 0) p[-del, ]
# 
# q = p
# q$tag = q$near
# q = split(q, list(q$clock, q$tag, q$intv))
# q = lapply(q, function(z) {
# 	if (nrow(z) == 0) return(NULL)
# 	del = which(z$conv < 1)
# 	if (length(del) > 0) z$conv[del] = 1e-8
# 	x = seq(log10(0.5), log10(20000), by=0.1)
# 	y = ecdf(log10(z$conv))(x)
# 	data.table(intv = z$intv[1], clock = z$clock[1], tag = z$tag[1], x=10^x, y=y)
# })
# q = rbindlist(q)
# 
# 
# conv = ggplot(q) +
# 	facet_grid(tag~clock) +
# 	geom_line(aes(x, y, color = factor(intv)), size = 1, alpha = 3/4) +
# 	scale_x_log10(breaks = 10^(0:5)) +
# 	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey60", "maroon1"))(3)[2],
# 																										 colorRampPalette(c("grey70", "navy"))(3)[2],
# 																										 colorRampPalette(c("grey50", "red"))(3)[2],
# 																										 colorRampPalette(c("grey95", "darkorange3"))(3)[2])))(length(unique(q$intv)))) +
# 	coord_cartesian(expand = F, xlim = c(0.5, 20000), ylim = c(-0.005, 1.005)) +
# 	theme_few() +
# 	theme(aspect.ratio = 1,
# 				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
# 				legend.title = element_text(size = 8.5),
# 				legend.text = element_text(size = 8),
# 				legend.key.height = unit(4, "mm"),
# 				legend.box.background = element_rect(fill = "white", size = 1/3, colour = "grey20"),
# 				legend.position = c(0.005, 0.99), legend.justification = c(c(0,1)), #legend.margin = margin(0,1,1,1, "mm"),
# 				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
# 				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
# 				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
# 				panel.grid.minor.x = element_blank(),
# 				panel.grid.major.y = element_line(size = 1/2, colour = "grey90"),
# 				panel.grid.minor.y = element_blank()) + 
# 	xlab("Width of 95% confidence interval (generations)") + ylab("CDF") +
# 	guides(color = guide_legend(title = "AF bin (%)", override.aes = list(size = 3, alpha = 1)))
# conv
# 
# ggsave(conv, filename = sprintf("_plot.conv.%s.pdf", make.names(tag)), height = 6.7, width = 10)

}





d = NULL

for (file in dir(pattern = "^dev_err_H_id_.+\\.age\\.pairs\\.txt$", path = ".", full.names = T)) {
	print(file)
	
	tag = sub("^dev_([^.]+)\\.age.+$", "\\1", basename(file))
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	tmp = cbind(tmp, times[tmp$MarkerID + 1, ])
	
	
	tmp$data = ""
	tmp$phased = ""
	tmp$expos = ""
	tmp$tmrca = as.numeric(sub("^.+_([0-9]+)_.+$", "\\1", tag))
	tmp$near = ""
	tmp$fltr = ""
	
	if (grepl("^tru_.+$", tag)) tmp$data = "(a) Before error"
	if (grepl("^err_.+$", tag)) tmp$data = "(b) After error"
	
	if (grepl("^.+_H_.+$", tag)) tmp$phased = "Correct phase"
	if (grepl("^.+_P_.+$", tag)) tmp$phased = "Phased"
	
	if (grepl("^.+_id_.+$", tag)) tmp$expos = "(a) Identical"
	if (grepl("^.+_fp_.+$", tag)) tmp$expos = "(b) False positives"
	if (grepl("^.+_fn_.+$", tag)) tmp$expos = "(c) False negatives"
	if (grepl("^.+_ff_.+$", tag)) tmp$expos = "(d) Both"
	
	if (grepl("_NN_", tag)) tmp$near = "Nearest neighbour"
	if (grepl("_RD_", tag)) tmp$near = "Random"
	
	if (grepl("_none$", tag)) tmp$fltr = "(a) No filter"
	if (grepl("_prop$", tag)) tmp$fltr = "(b) Fixed proportion"
	if (grepl("_auto$", tag)) tmp$fltr = "(c) Minimum exclusion"
	
	tmp$Clock[which(tmp$Clock == "M")] = "(a) Mutation clock"
	tmp$Clock[which(tmp$Clock == "R")] = "(b) Recombination clock"
	tmp$Clock[which(tmp$Clock == "C")] = "(c) Combined clock"
	
	tmp$tag = tag
	
	
	tmp$i.len = (tmp$SegmentRHS - tmp$SegmentLHS) + 1
	tmp$p.len = (marker$Position[tmp$SegmentRHS+1] - marker$Position[tmp$SegmentLHS+1]) + 1
	tmp$g.len = (marker$GenDist[tmp$SegmentRHS+1] - marker$GenDist[tmp$SegmentLHS+1]) + 1e-8
	
	
	d = rbind(d, tmp)
}


#d = d[(which(d$phased == "Correct phase")), ]


d$p.len = d$p.len / 1e6

d$cordant = "Discordant"
d$cordant[which(d$Shared == 1)] = "Concordant"


#del = which(d$SegmentLHS <= 1 | d$SegmentRHS >= max(d$SegmentRHS)-1)
#if (length(del) > 0) d = d[-del,]



#d$expos = factor(d$expos, levels = sort(unique(d$expos)), labels = sub("^\\([a-z]\\) (.+)", "\\1", sort(unique(d$expos))))


#D = split(d, list(d$data, d$expos))
p = d

#for (tag in names(D)) {
#print(tag)

#p = D[[tag]]

nbins = 11

p$intv = ((cut(p$Fk / 5000 * 100, breaks = (0:10)/10, include.lowest = T)))
if (any(is.na(p$intv))) {
	p$intv = ((cut(p$Fk / 5000 * 100, breaks = (0:11)/10, include.lowest = T)))
	nbins = 12
}

gr = NULL
for (k in 2:52) gr = rbind(gr, data.table(x=k-0.5, y = 10^(-8:5)))

xl = (1:(length(levels(p$intv))))/10
xl[ 1:length(xl) %% 2 == 0 ] = ""


len = ggplot(p) +
	facet_grid(near+fltr~Clock+cordant) +
	#geom_raster(aes(intv,    g.len), position = "identity", stat = "bin2d", binwidth = c(1/11,1/5)) +
	geom_tile(aes(intv,      g.len, fill=..density..*100), position = "identity", stat = "bin2d", binwidth = c(1/nbins,1/5)) +
	#geom_tile(data = q, aes(intv, ymin=ymin, ymax=ymax, fill=value), position = "identity") +
	geom_vline(xintercept = (1:(nbins-1))-0.5, alpha = 1/10, color = "black") +
	#geom_point(data = gr, aes(x, y), alpha = 1/2, color = "black", shape = 95, size = 2) +
	geom_linerange(aes(intv, g.len), stat = "summary", alpha = 1/4, color = "black", fun.ymin = function(z) {min(z)}, fun.ymax = function(z) {max(z)}) +
	geom_linerange(aes(intv, g.len), stat = "summary", size = 1.5, alpha = 1/2, color = "black", fun.ymin = function(z) {quantile(z,0.25)}, fun.ymax = function(z) {quantile(z,0.75)}) +
	geom_point(aes(intv,     g.len), stat = "summary", alpha = 2/4, color = "black", fun.y = median, shape=19) +
	geom_point(aes(intv,     g.len), stat = "summary", alpha = 2/4, color = "white", fun.y = median, shape=20) +
	coord_cartesian(ylim = c(7e-09, 150), expand = F) +
	#scale_x_continuous(breaks = (1:10)/10) +
	scale_x_discrete(labels = xl) +
	scale_y_log10(breaks = 10^(-8:2), labels = c("", trimws(format(10^(-7:1), scientific = T)), "")) +
	scale_fill_gradientn(colours = c("white", "goldenrod1", "orange", "red", "darkred"), na.value = "darkred", limits=c(0, 25), breaks = c(0,5,10,15,20,25), labels=sprintf("%d%%", c(0,5,10,15,20,25))) +
	theme_few() +
	theme(aspect.ratio = 2,
				axis.text = element_text(size = 8),
				axis.text.x = element_text(size = 6),
				panel.border=element_rect(fill = NA, colour = "grey50", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,0), legend.position=c(1-0.995,1-0.995),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2),
				#legend.position = "bottom",
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 7),
				legend.title = element_blank()) +
	ylab("Genetic length (cM)") + 
	xlab("Allele frequency bin (%)")
#len

ggsave(len, filename = sprintf("_plot_pairs.len-g.%s.pdf", make.names(tag)), height = 12.5, width = 7.5)


len = ggplot(p) +
	facet_grid(near+fltr~Clock+cordant) +
	#geom_raster(aes(intv,    g.len), position = "identity", stat = "bin2d", binwidth = c(1/11,1/5)) +
	geom_tile(aes(intv,      p.len, fill=..density..*100), position = "identity", stat = "bin2d", binwidth = c(1/nbins,1/5)) +
	#geom_tile(data = q, aes(intv, ymin=ymin, ymax=ymax, fill=value), position = "identity") +
	geom_vline(xintercept = (1:(nbins-1))-0.5, alpha = 1/10, color = "black") +
	#geom_point(data = gr, aes(x, y), alpha = 1/2, color = "black", shape = 95, size = 2) +
	geom_linerange(aes(intv, p.len), stat = "summary", alpha = 1/4, color = "black", fun.ymin = function(z) {min(z)}, fun.ymax = function(z) {max(z)}) +
	geom_linerange(aes(intv, p.len), stat = "summary", size = 1.5, alpha = 1/2, color = "black", fun.ymin = function(z) {quantile(z,0.25)}, fun.ymax = function(z) {quantile(z,0.75)}) +
	geom_point(aes(intv,     p.len), stat = "summary", alpha = 2/4, color = "black", fun.y = median, shape=19) +
	geom_point(aes(intv,     p.len), stat = "summary", alpha = 2/4, color = "white", fun.y = median, shape=20) +
	coord_cartesian(ylim = c(7e-09, 150), expand = F) +
	#scale_x_continuous(breaks = (1:10)/10) +
	scale_x_discrete(labels = xl) +
	scale_y_log10(breaks = 10^(-8:2), labels = c("", trimws(format(10^(-7:1), scientific = T)), "")) +
	scale_fill_gradientn(colours = c("white", "goldenrod1", "orange", "red", "darkred"), na.value = "darkred", limits=c(0, 25), breaks = c(0,5,10,15,20,25), labels=sprintf("%d%%", c(0,5,10,15,20,25))) +
	theme_few() +
	theme(aspect.ratio = 2,
				axis.text = element_text(size = 8),
				axis.text.x = element_text(size = 6),
				panel.border=element_rect(fill = NA, colour = "grey50", size = 2/3),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,0), legend.position=c(1-0.995,1-0.995),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2),
				#legend.position = "bottom",
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 7),
				legend.title = element_blank()) +
	ylab("Physical length (Mb)") + 
	xlab("Allele frequency bin (%)")
#len

ggsave(len, filename = sprintf("_plot_pairs.len-p.%s.pdf", make.names(tag)), height = 12.5, width = 7.5)







q = data.table(MarkerID = p$MarkerID, Clock = p$Clock, y = p$q50, tc = p$node.time, td = p$parent.time, cordant = p$cordant, 
							 data = p$data, phased = p$phased, expos = p$expos, tmrca = p$tmrca, near = p$near, fltr = p$fltr, tag = p$tag)

q = split(q, sprintf("%d %s %s %s", q$MarkerID, q$Clock, q$cordant, q$tag))
q = lapply(q, function(x) {
	if (x$cordant[1] == "Concordant") {
		i = which.max(x$y)
		x = x[i,]
	} else {
		i = which.min(x$y)
		x = x[i,]
	}
	x
})
q = rbindlist(q)

w = split(q, sprintf("%s %s %s", q$Clock, q$cordant, q$tag))
w = lapply(w, function(x) {
	if (x$cordant[1] == "Concordant") {
		x$rankc = cor(x$tc, x$y, method = "spearman")
		x$rmsle = rmsle(x$tc + 1e-32, x$y * Ne + 1e-32)
	} else {
		x$rankc = cor(x$td, x$y, method = "spearman")
		x$rmsle = rmsle(x$td + 1e-32, x$y * Ne + 1e-32)
	}
	x[1,]
})
w = rbindlist(w)

w$rankc.txt = sprintf("Spearman's rho = %.3f", w$rankc)
w$rmsle.txt = sprintf("RMSLE = %.3f", w$rmsle)


q$x = q$tc
x = which(q$cordant == "Discordant")
q$x[x] = q$td[x]


mx = c(0.5, 15e4)

tm = rbind(data.table(x = 10^(0:6), y = 1e-6, xend = 10^(0:6), yend = 0.7),
					 data.table(y = 10^(0:6), x = 1e-6, yend = 10^(0:6), xend = 0.7))

tn = c((2:9), 10*(2:9), 100*(2:9), 1000*(2:9), 10000*(2:9))
tn = rbind(data.table(x = tn, xend = tn, y = 1e-6, yend = 0.615), 
					 data.table(y = tn, yend = tn, x = 1e-6, xend = 0.615))

tx = rbind(data.table(x = 1e-6, xend = 1e6, y = 1e-6, yend = 0.7),
					 data.table(x = 1e-6, xend = 0.7, y = 1e-6, yend = 1e6))


scat = ggplot(q) +
	facet_grid(near+fltr~cordant+Clock) +
	geom_raster(aes(x, y * Ne), position = "identity", stat = "bin2d", binwidth = c(0.075, 0.075)) +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_smooth(aes(x, y * Ne), method = "lm", color = "black", size = 1/2, se = F, fullrange = F) +
	geom_smooth(aes(x, y * Ne), method = "lm", color = "white", linetype = "22", size = 1/2, fill = "black", se = T, fullrange = F) +
	geom_text(data = w, aes(x = 40, y = 80000, label = rankc.txt), size = 2.5) +
	geom_text(data = w, aes(x = 8000, y = 1.5, label = rmsle.txt), size = 2.5) +
	geom_rect(data = tx, aes(xmin=x, xmax=xend, ymin=y, ymax=yend), fill="grey80") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), size = 1/2, color="black") +
	geom_segment(data = tn, aes(x=x, xend=xend, y=y, yend=yend), size = 1/4, color="black") +
	coord_cartesian(xlim = mx, ylim = mx, expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("lightcyan1", "lightblue1", "lightblue2", "deepskyblue", "dodgerblue", "dodgerblue2", "purple2", "purple3", "purple4"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey90", "grey80", "lightskyblue", "deepskyblue", "dodgerblue", "purple1", "purple4", "navy", "black"), na.value = "black", trans = "log10", limits = c(1, 500)) +
	#scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 6),
				#axis.ticks = element_blank(),
				panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				legend.justification=c(0,1), legend.position=c(1-0.99,0.975),
				legend.background = element_rect(fill = NA, colour = NA, size = NA),
				#legend.position = "bottom",
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 6),
				legend.title = element_blank()) +
	ylab("Inferred TMRCA") + xlab("True TMRCA")
#scat

ggsave(scat, filename = sprintf("_plot_pairs.scat.%s.pdf", make.names(tag)), height = 12.5, width = 12.5)

	



tmp = split(p, list(p$Clock, p$intv, p$fltr, p$near))
tmp = lapply(tmp, function(x) {
	if (nrow(x) == 0) return(NULL)
	
	con = which(x$Shared == 1)
	dis = which(x$Shared == 0)
	
	ic = seq(-1.75, 0.25, by=0.01)
	id = seq(-0.25, 1.75, by=0.01)
	
	mc = log10(((x$q50[con] * Ne) / (x$node.time[con])))
	mc = ecdf(mc)(ic)
	
	md = log10(((x$q50[dis] * Ne) / (x$parent.time[dis])))
	md = ecdf(md)(id)
	
	rbind(data.table(intv = x$intv[1], clock = x$Clock[1], tag1 = x$fltr[1], tag2 = x$near[1], type = "Concordant", x=ic, y=mc),
				data.table(intv = x$intv[1], clock = x$Clock[1], tag1 = x$fltr[1], tag2 = x$near[1], type = "Discordant", x=id, y=md))
})
q = rbindlist(tmp)


hist = ggplot(q) +
	facet_grid(tag2+tag1~type+clock, scales = "free_x") +
	geom_vline(xintercept = 0, linetype = "dashed") +
	geom_line(aes(x, y, color = (intv))) +
	scale_y_continuous(breaks = (1:9)/10) +
	scale_color_manual(values = colorRampPalette(rev(c(colorRampPalette(c("grey60", "maroon1"))(3)[2],
																										 colorRampPalette(c("grey70", "navy"))(3)[2],
																										 colorRampPalette(c("grey50", "red"))(3)[2],
																										 colorRampPalette(c("grey95", "darkorange3"))(3)[2])))(length(unique(q$intv)))) +
	coord_cartesian(expand = F, ylim = c(-0.005, 1.005)) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title = element_text(size = 7.5),
				legend.text = element_text(size = 7),
				legend.key.height = unit(3.5, "mm"),
				legend.key.width = unit(3.5, "mm"),
				legend.box.background = element_rect(fill = "white", size = 1/4, colour = "grey40"),
				legend.position = c(0.0025, 0.9975), legend.justification = c(c(0,1)), legend.margin = margin(1,1,1,1, "mm"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.y = element_blank()) + 
	xlab("Relative pairwise TMRCA estimate, median (log-scale)") + ylab("CDF") +
	guides(color = guide_legend(title = "AF bin (%)", override.aes = list(size = 2, alpha = 1)))
#hist

ggsave(hist, filename = sprintf("_plot_pairs.hist.%s.pdf", make.names(tag)), height = 12.5, width = 12.5)





q = data.table(MarkerID = p$MarkerID, Clock = p$Clock, y = p$q50, tc = p$node.time, td = p$parent.time, cordant = p$cordant, intv = p$intv,
							 data = p$data, phased = p$phased, expos = p$expos, tmrca = p$tmrca, near = p$near, fltr = p$fltr, tag = p$tag)

q = split(q, sprintf("%d %s %s", q$MarkerID, q$Clock, q$tag))
q = lapply(q, function(x) {
	if (nrow(x) == 0) return(NULL)
	
	con = which(x$cordant == "Concordant")
	dis = which(x$cordant == "Discordant")
	
	if (length(con) == 0) return(NULL)
	if (length(dis) == 0) return(NULL)
	
	con = x[con, ]
	dis = x[dis, ]
	
	i = which.max(con$y)
	con = con[i, ]
	
	i = which.min(dis$y)
	dis = dis[i, ]
	
	rbind(con, dis)
})
q = rbindlist(q)

w = split(q, sprintf("%s %s %s %s", q$Clock, q$intv, q$cordant, q$tag))
w = lapply(w, function(x) {
	a = x
	b = x
	a$method = "rankc"
	b$method = "rmsle"
	if (x$cordant[1] == "Concordant") {
		a$val = cor(x$tc, x$y, method = "spearman")
		b$val = rmsle(x$tc + 1e-32, x$y * Ne + 1e-32)
	} else {
		a$val = cor(x$td, x$y, method = "spearman")
		b$val = rmsle(x$td + 1e-32, x$y * Ne + 1e-32)
	}
	rbind(a[1, ], b[1, ])
})
w = rbindlist(w)

xl = (1:(length(levels(w$intv))))/10

w = split(w, w$method)



aa = ggplot(w$rankc) +
	facet_grid(cordant~Clock, drop = T) +
	geom_line(aes(as.numeric(intv), val, color = factor(fltr), linetype = factor(near)), alpha = 0.8) +
	geom_point(aes(intv, val, color = factor(fltr), shape = factor(near)), alpha = 0.8) +
	coord_cartesian(ylim=c(-0.01, 1.01), xlim = c(0.5, 10.5), expand = F) +
	scale_x_discrete(labels = xl) +
	scale_y_continuous(breaks = (1:20)/5) +
	theme_bw() +
	theme(#aspect.ratio = 2/3,
		axis.text = element_text(size = 8),
		panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
		strip.text.x = element_text(face = "bold", size = 9),
		strip.text.y = element_text(face = "bold", angle = -90, size = 9),
		#legend.justification=c(1,1), legend.position=c(0.99,0.99),
		#legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
		legend.position = "bottom",
		legend.text = element_text(size = 9),
		legend.title = element_blank()) +
	xlab("Allele frequency bin (%)") + ylab("Spearman's rho")
aa

ggsave(aa, filename = sprintf("_plot_pairs.stat_rankc.%s.pdf", make.names(tag)), height = 7.5, width = 7.5)


bb = ggplot(w$rmsle) +
	facet_grid(cordant~Clock, drop = T) +
	geom_line(aes(as.numeric(intv), val, color = factor(fltr), linetype = factor(near)), alpha = 0.8) +
	geom_point(aes(intv, val, color = factor(fltr), shape = factor(near)), alpha = 0.8) +
	coord_cartesian(ylim=c(-0.01, 1.95), xlim = c(0.5, 10.5), expand = F) +
	scale_x_discrete(labels = xl) +
	scale_y_continuous(breaks = (1:20)/5) +
	theme_bw() +
	theme(#aspect.ratio = 2/3,
		axis.text = element_text(size = 8),
		panel.border = element_rect(fill = NA, colour = "grey50", size = 1/2),
		strip.text.x = element_text(face = "bold", size = 9),
		strip.text.y = element_text(face = "bold", angle = -90, size = 9),
		#legend.justification=c(1,1), legend.position=c(0.99,0.99),
		#legend.background = element_rect(fill = "white", colour = "grey80", size = 1/3),
		legend.position = "bottom",
		legend.text = element_text(size = 9),
		legend.title = element_blank()) +
	xlab("Allele frequency bin (%)") + ylab("RMSLE")
bb

ggsave(bb, filename = sprintf("_plot_pairs.stat_rmsle.%s.pdf", make.names(tag)), height = 7.5, width = 7.5)






}










a = sapply(split(p, list(p$clock)), function(x) {
	data.table("tc" = round(cor(x$PostMode * Ne, x$node.time, method = "s"), 4),
						 "tm" = round(cor(x$PostMode * Ne, x$mid, method = "s"), 4),
						 "td" = round(cor(x$PostMode * Ne, x$parent.time, method = "s"), 4))
	
})
a = t(array(a, dim = dim(a), dimnames = dimnames(a)))
a

b = sapply(split(p, list(p$clock)), function(x) {
	data.table("tc" = round(rmsle(x$PostMode * Ne, x$node.time), 4),
						 "tm" = round(rmsle(x$PostMode * Ne, x$mid), 4),
						 "td" = round(rmsle(x$PostMode * Ne, x$parent.time), 4))
	
})
b = t(array(b, dim = dim(b), dimnames = dimnames(b)))
b


x = by(p, list(cut(p$Fk, breaks = c(2, 5, 10, 20, 30, 40, 50), include.lowest = T), p$clock), function(x) cor(x$PostMode * Ne, x$parent.node, method = "s"))
x = melt(array(x, dim = dim(x), dimnames = dimnames(x)))

ggplot(x) +
	geom_bar(aes(x=Var1, y=value, fill = Var2), stat = "identity", position = "dodge")

x = by(p, list(p$Fk, p$clock, p$phased), function(x) rmsle(x$PostMode * Ne, x$node.time))
array(x, dim = dim(x), dimnames = dimnames(x))


x = by(p, list(p$Fk, p$clock, p$phased), nrow)
array(x, dim = dim(x), dimnames = dimnames(x))













q = data.table(clock = d$clock,
							 method = d$method,
							 con = d$Lower * Ne,
							 dis = d$Upper * Ne,
							 est = d$PostMode * Ne,
							 t0 = d$node.time,
							 t1 = d$parent.time,
							 err = d$error)

q$est = log(q$est / q$t0) / log(q$t1 / q$t0)


ggplot(q) +
	facet_grid(method~clock) +
	#geom_rect(data = tr, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), size = 1, color = "black", fill = NA) +
	# geom_vline(xintercept = 0, size = 3/4, linetype = "21", color = "green4") +
	# geom_vline(xintercept = 1, size = 3/4, linetype = "21", color = "sienna") +
	geom_vline(xintercept = c(0,1), size = 3/4, linetype = "21") +
	geom_vline(xintercept = 0.5, size = 3/4, color = "grey60") +
	geom_histogram(data = q[which(q$err == F), ], aes(est), binwidth = 0.1, alpha = 1/2, fill = "seagreen") +
	geom_histogram(data = q[which(q$err == T), ], aes(est), binwidth = 0.1, alpha = 1/2, fill = "red3") +
	scale_x_continuous(breaks = -10:10) +
	scale_y_continuous(breaks = seq(0, 10000, by = 250), labels = format(seq(0, 10000, by = 250), big.mark = ",")) +
	coord_cartesian(expand = F, xlim = c(-3.5, 4.5), ylim = c(0, 1550)) +
	theme_few() +
	theme(aspect.ratio = 1,
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 9),
				strip.text.y = element_text(face = "bold", angle = -90, size = 11),
				panel.grid.major.x = element_line(size = 1/2, colour = "grey90"),
				panel.grid.minor.x = element_line(size = 1/2, colour = "grey90")) +
	xlab("Relative age estimate (log-scale)") + ylab("Count")

ggsave(filename = "_plot.hist.pdf", height = 10, width = 8, dpi = 300)


z = by(p, list(p$clock, p$method), function(p) length(which(p$PostMode*Ne < p$node.time)) / nrow(p) * 100)
array(z, dim = dim(z), dimnames = dimnames(z))
z = by(p, list(p$clock, p$method), function(p) length(which(p$PostMode*Ne > p$parent.time)) / nrow(p) * 100)
array(z, dim = dim(z), dimnames = dimnames(z))




q = rbindlist(lapply(split(q, q$clock), function(p) {
	tmp = NULL
	for (tag in names(p)[-1]) {
		den = density(p[[tag]], from = -5, to = 6, n = 100)
		tmp = rbind(tmp, data.table(clock = p$clock[1], type = tag, x = den$x, y = den$y))
	}
	tmp
}))


ggplot(q) +
	facet_wrap(~clock, ncol = 1) +
	geom_vline(xintercept = c(0, 1)) +
	#geom_area(aes(x=x, y = y, fill = type), alpha = 1/3, position = position_identity()) +
	#geom_line(aes(x=x, y = y, color = type), alpha = 2/3, size = 1) +
	geom_bar(aes(x=x, y=y), stat = "identity") +
	#geom_hline(yintercept = 0, color = "white", size = 1) +
	#scale_colour_manual(values = c(con = "forestgreen", dis = "saddlebrown", est = "dodgerblue")) +
	#scale_fill_manual(values = c(con = "forestgreen", dis = "saddlebrown", est = "dodgerblue1")) +
	scale_x_continuous(breaks = -10:10) +
	coord_cartesian(expand = F, xlim = c(-3.5, 4.5), ylim = c(0, max(q$y + 1))) +
	theme_few() +
	theme(aspect.ratio = 1/2)



ggplot(p) +
	#facet_grid(fk ~ .) +
	#geom_histogram(aes(est), binwidth = 0.05, alpha = 0.5, fill = "deepskyblue") +
	geom_density(aes(con, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "forestgreen", fill = "green4") +
	stat_bin(aes(dis, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "saddlebrown") +
	stat_bin(aes(est, ..density..), geom = "step", binwidth = 0.05, size = 1, alpha = 0.75, color = "deepskyblue") +
	geom_vline(xintercept = c(0, 1)) +
	scale_x_continuous(breaks = -10:10) +
	coord_cartesian(xlim = c(-3, 4)) +
	theme_few()



#ggplot(p) + geom_point(aes(time, PostMode * Ne))


ggplot(p) +
	geom_smooth(aes(time, PostMode * Ne), method = "lm") +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_linerange(aes(x = time, ymax = Upper * Ne, ymin = Lower * Ne)) +
	geom_point(aes(time, PostMode * Ne)) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F)))


#d$time = times$time[ d$MarkerID + 1 ]


d$method[which(d$method == "FGT")] = "bold('(a) FGT, true haplotypes')"
d$method[which(d$method == "DGT")] = "bold('(c) DGT')"
d$method[which(d$method == "HMM")] = "bold('(d) HMM')"

d$method[which(d$type == "truP" | d$type == "errP")] = "bold('(b) FGT, phased haplotypes')"

d$clock[which(d$clock == "Mut")]    = " bolditalic(T[M]) "
d$clock[which(d$clock == "Rec")]    = " bolditalic(T[R]) "
d$clock[which(d$clock == "MutRec")] = "bolditalic(T[MR])"

d$type[which(d$type == "truH" | d$type == "truP")] = "tru"
d$type[which(d$type == "errH" | d$type == "errP")] = "err"


d$type = sprintf("bold('%s')", d$type)


p = d[which(d$type == "tru"), ]
p = d[which(d$type == "err"), ]

# subset
q = split(p, list(p$method, p$clock))
q = lapply(q, function(p) {
	if (nrow(p) < 5000) return(p)
	p[sample(1:nrow(p), 5000), ]
})
p = rbindlist(q)


tm = rbind(data.table(x = 10^(0:6), y = 0.01, xend = 10^(0:6), yend = 0.65),
					 data.table(y = 10^(0:6), x = 0.01, yend = 10^(0:6), xend = 0.65))

gg = ggplot(p) +
	facet_grid(clock~type, labeller = label_parsed) +
	# geom_rect(data = data.table(xmin = 0.2, xmax = 50000, ymin = 0.2, ymax = 50000),
	# 					aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey80") +
	geom_raster(aes(time, PostMedian * Ne), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_bin2d(aes(time, PostMean * 20000), binwidth = c(0.075, 0.075), size = 2) +
	#stat_density2d(aes(time, PostMean * 10000, fill = ..level.., alpha = ..level..), size = 0.01, bins = 50, geom = 'polygon') +
	geom_density_2d(aes(time, PostMedian * Ne), bins = 5, size = 1/3, color = "black", alpha = 0.8) +
	geom_smooth(aes(time, PostMedian * Ne), method = "lm") +
	geom_abline(intercept = c(0,0), slope = 1, alpha = 1/3) +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	coord_cartesian(xlim = c(0.5, 50000), ylim = c(0.5, 50000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("white", "deepskyblue", "royalblue2", "magenta4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("grey75", "yellow", "orange", "orangered", "darkred"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	#scale_color_gradientn(colours = c("lightblue", "deepskyblue", "royalblue", "darkblue"), na.value = "navy", limits = c(0, 4)) +
	#scale_fill_gradient(low = "green", high = "red") +
	scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				#panel.spacing.y = unit(0.05, "cm"),
				#plot.margin = unit(c(-0.1, 0.05, 0.05, 0.05), "cm"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				#panel.background = element_rect(fill = "grey80"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 8),
				strip.text.y = element_text(face = "bold.italic", angle = 0),
				legend.justification=c(1,0), legend.position=c(0.995,1-0.995), legend.box.background = element_rect(colour = "black"),
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age") #+ labs(fill = "Density")
gg

ggsave(gg, filename = "_plot.age.check.HMM.postmean.pdf", width = 6, height = 8)
ggsave(gg, filename = "_plot.age.check.HMM.postmode.pdf", width = 6, height = 8)
ggsave(gg, filename = "_plot.age.check.HMM.postmedi.pdf", width = 6, height = 8)

ggsave(gg, filename = "_plot.age.sites.africa.tru.pdf", width = 8, height = 6.5)
ggsave(gg, filename = "_plot.age.sites.africa.err.pdf", width = 8, height = 6.5)



q = split(d, list(d$type, d$method, d$clock, d$ndiscord))
q = lapply(q, function(p) {
	if (nrow(p) < 5000) return(p)
	p[sample(1:nrow(p), 5000), ]
})
d = rbindlist(q)


x = by(d, list(d$Fk, d$type, d$method, d$clock), function(x) cor(x$time, x$PostMean, method = "s"))
array(x, dim = dim(x), dimnames = dimnames(x))

y = by(d, list(d$Fk, d$type, d$method, d$clock), function(x) rmsle(x$time, x$PostMedian * Ne))
array(y, dim = dim(y), dimnames = dimnames(y))

z = by(d, list(d$type, d$method, d$clock), function(x) mean(abs(x$time - x$PostMean * Ne)))
array(z, dim = dim(z), dimnames = dimnames(z))





x = by(d, list(d$type, d$method, d$clock), nrow)
array(x, dim = dim(x), dimnames = dimnames(x))



ggplot(p) +
	#facet_grid(clock~., labeller = label_parsed) +
	stat_summary(aes(Fk, PostMean * 20000), fun.data = "median_hilow", geom = "smooth", alpha = 0.5) +
	stat_summary(aes(Fk, time), fun.data = "median_hilow", geom = "smooth", color = "black", alpha = 0.5) +
	scale_y_log10() +
	scale_x_log10() +
	theme_few()




###
### ration of con/dis lengths
###

load("data.result.age.pairs.check.RData")

markers = fread("vanilla.marker.txt", header = T)

d$len = markers$Position[d$SegmentRHS + 1] - markers$Position[d$SegmentLHS + 1] + 1

d$Shared = as.character(d$Shared)

ggplot(d) +
	#facet_grid(clock~type) +
	geom_histogram(aes(len, fill = Shared), bins = 100, position=position_dodge()) +
	scale_x_log10() +
	scale_y_log10() +
	theme_bw() +
	xlab("Physical length") + ylab("Count")




###












markers = fread("vanilla.marker.txt", header = T, stringsAsFactors = F)

files = dir(pattern = ".+\\.age\\.sites\\..+")



file = "test.age.sites.FGT.MutRec.txt"
file = "test.age.sites.FGT.Mut.txt"
file = "test.age.sites.FGT.Rec.txt"

file = "test.age.sites.DGT.MutRec.txt"
file = "test.age.sites.DGT.Mut.txt"
file = "test.age.sites.DGT.Rec.txt"

file = "test.age.sites.HMM.MutRec.txt"
file = "test.age.sites.HMM.Mut.txt"
file = "test.age.sites.HMM.Rec.txt"

age = read.age(file, M, R)

ggplot(age) +
	geom_point(aes(time, PostMean * 10000, color = (Fk))) +
	geom_abline(intercept = c(0,0), slope = 1) +
	coord_cartesian(xlim = c(0.2, 200000), ylim = c(0.2, 200000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = format(10^(0:6), big.mark = ',', scientific = F)) +
	scale_y_log10(breaks = 10^(0:6), labels = format(10^(0:6), big.mark = ',', scientific = F))

cor(age$time, age$PostMean, method = "s")
cor(age$Fk, age$PostMean, method = "s")
cor(age$time, age$Fk, method = "s")

cor(age$time, age$PostMedian, method = "s")
cor(age$time, age$PostMode, method = "s")
cor(age$time, age$Robust, method = "s")
cor(age$time, age$Lower, method = "s")

plot(sort(age$PostMean*20000), ylim=c(0, 20000))
plot(sort(age$time), ylim=c(0,20000))

cor(age$PostMean, age$Fk)^2

plot(age$Fk, age$PostMean)


rmsle = function(a, b) {
	n = length(a)
	a = log(a + 1)
	b = log(b + 1)
	sqrt(sum((a - b)^2) / n)
}

rmsle(age$time, age$PostMean)


d = NULL

for (file in files) {
	
	method = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\1", file)
	clock  = sub(".+\\.age\\.sites\\.(.+)\\.(.+)\\.txt", "\\2", file)
	
	phased = grepl(".+P\\.age\\..+", file)
	
	if (phased)
		method = sprintf("%s phased", method)
	
	#if (clock == "MutRec") clock = " MutRec"
	
	cat(".")
	
	age = read.age(file, M, R)
	
	age$method = method
	age$clock = clock
	
	d = rbind(d, age)
}


del = which(d$N_Others < d$N_Shared)
if (length(del) > 0)
	d = d[-del, ]


m = match(d$MarkerID, markers$MarkerID)
d$a = markers$AlleleCount1[m]
d$g = markers$GenotypeCount1[m]
del = which(d$a != d$g)
if (length(del) > 0)
	d = d[-del, ]


#d = d[sample(1:nrow(d), 50000), ]


d$frq = (d$Fk / 5000) * 100

d$est = NA

x = (d$clock == "MutRec")
d$est[x] = d$PostMean[x] * 7300 * (2/100) * 2
x = (d$clock == "Mut")
d$est[x] = d$PostMean[x] * 7300 * (2/100)
x = (d$clock == "Rec")
d$est[x] = d$PostMean[x] * 7300 * 2

x = (d$method == "FGT" | d$method == "FGT phased")
d$est[x] = d$PostMean[x] * 7300 * (1/4)


x = by(d, list(d$clock, d$method), function(x) median(abs(x$time - x$est)) )
array(x, dim(x), dimnames(x))

x = by(d, list(d$clock, d$method), function(x) cor(x$time, x$est, method = "s") )
array(x, dim(x), dimnames(x))



gg = ggplot(d) +
	facet_grid(clock~method) +
	geom_point(aes(x=frq, y=est), colour = "grey50", alpha = 0.05, shape=15, size=1) +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=2.5, colour = "white") +
	stat_summary(aes(x=frq, y=est), geom = "line", fun.y = "median", size=1, colour = "grey10") +
	stat_summary(aes(x=frq, y=time), geom = "point", fun.y = "median", size=1.5, colour = "blue") +
	coord_cartesian(ylim = c(1, 2000)) +
	scale_y_log10(breaks=c(0, 1, 10, 100, 1000, 10000, 100000)) +
	scale_x_continuous(breaks = (0:5)/10) +
	theme_few() +
	theme(aspect.ratio=2.1/3) +
	xlab("Allele frequency (%)") + ylab("Inferred age (generations)")

ggsave(gg, filename = "_plot.ooa.tru.png", width = 12, height = 10)
ggsave(gg, filename = "_plot.ooa.err.png", width = 12, height = 10)



ggplot(d) +
	facet_grid(clock~method) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(time, est, colour = Fk), size = 0.5) +
	scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	scale_y_log10(breaks=c(1, 10, 100, 1000, 10000, 100000)) +
	coord_cartesian(xlim = c(0.05, 50000), ylim = c(0.05, 50000)) +
	theme(aspect.ratio = 1)  +
	xlab("True age") + ylab("Posterior mode age")





files = dir(pattern = "^t1000\\.truH\\.age\\.pairs\\..+")

p = lapply(files, function(file) {
	x = fread(file, header = T, stringsAsFactors = F)
	x$method = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\1", file)
	x$clock  = sub(".+\\.age\\.pairs\\.(.+)\\.(.+)\\.txt", "\\2", file)
	x
})
p = rbindlist(p)

p = p[which(p$Shared == 1),]
p = p[-(which(p$SegmentLHS == 0)),]
p = p[-(which(p$SegmentRHS == max(p$SegmentRHS))),]

pos = markers$Position
names(pos) = as.character(markers$MarkerID)

p$pos_SegmentLHS = pos[ as.character(p$SegmentLHS) ]
p$pos_SegmentRHS = pos[ as.character(p$SegmentRHS) ]

p$len = (p$pos_SegmentRHS - p$pos_SegmentLHS) + 1

p$str = sprintf("%d %d %d", p$MarkerID + 1, p$SampleID0 + 1, p$SampleID1 + 1)



load("~/Research/DetectIBD/result.truth.local.RData")

truth$str = sprintf("%d %d %d", truth$index, truth$g0, truth$g1)


p = split(p, list(p$method, p$clock))

for (tag in names(p)) {
	sub = p[[tag]]
	
	i = intersect(sub$str, truth$str)
	
	a = match(i, sub$str)
	b = match(i, truth$str)
	
	sub = sub[a, ]
	t = truth[b, ]
	
	sub$len_true = (t$rhs.position - t$lhs.position) + 1
	
	p[[tag]] = sub
}

p = rbindlist(p)



ggplot(p[sample(1:nrow(p), 100000),]) +
	facet_grid(method~clock) +
	geom_point(aes(len_true, len)) +
	scale_x_log10() +
	scale_y_log10()
stat_summary(aes(x=Fk, y=time), geom = "line", fun.y = "median", size=2, colour="white", alpha = 0.7) +
	geom_boxplot(aes(factor(Fk), len))








M = read.table("vanilla.mutations.txt", header=T)
R = read.table("vanilla.records.txt", header=T)

get.age = function(age, M, R) {
	age$node = M$node[ age$MarkerID + 1 ]
	
	x = match(age$node, R$node)
	
	del = which(is.na(x))
	if (length(del) > 0)
		age = age[-del, ]
	
	age$time = R$time[ x ]
	age$length = R$right[ x ] - R$left[ x ]
	
	return(age)
}

read.age = function(file, M, R) {
	age = fread(file, header = T, stringsAsFactors = F)
	get.age(age, M, R)
}






z = fread("x2349.100.age.pairs.FGT.MutRec.HardBreaks.txt", header = T)

z$node = M$node[ z$MarkerID + 1 ]

x = match(z$node[1], R$node)


c0 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\1", R$children))
c1 = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\2", R$children))

a = match(60*2+2, c0)
b = match(75*2, c1)
intersect(a, b)


p = data.table(fk  = d$Fk,
							 shr = d$Shared,
							 lhs = M$position[z$SegmentLHS+1],
							 rhs = M$position[z$SegmentRHS+1])
p$len = p$rhs - p$lhs + 1

save(d, file="data.age-ibd.pairs.RData")

ggplot(d) + facet_wrap(~shr) + stat_summary(aes(fk, len), fun.data = "median_hilow", geom = "smooth") + scale_y_log10()


###
stop()
###

library(data.table)

files = dir(pattern="^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$ndiscord = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	tmp$method = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\3", file)
	tmp$clock = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.age.sites.HardBreaks.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$ndiscord = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	tmp$method = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\3", file)
	tmp$clock = sub("^y([0-9]+)\\.([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.age.sites.SoftBreaks.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^x([0-9]+)\\.1000\\.age\\.pairs\\.([A-Z]+)\\.MutRec\\.HardBreaks\\.txt$")

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$method = sub("^x([0-9]+)\\.([0-9]+)\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\3", file)
	d = rbind(d, tmp)
}

save(d, file="data.age-ibd.pairs.RData")


M = read.table("vanilla.mutations.txt", header=T)

p = data.table(fk  = d$Fk,
							 shr = d$Shared,
							 lhs = M$position[d$SegmentLHS+1],
							 rhs = M$position[d$SegmentRHS+1])
p$len = p$rhs - p$lhs + 1

save(p, file="data.age-ibd.pairs.parsed.RData")




###
stop()
###

library(data.table)

files = dir(pattern="^result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", recursive = T, full.names = T)

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	tmp$type   = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
	tmp$method = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\3", file)
	tmp$clock  = sub("^./(.+)/result\\.sample([0-9]+)\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\4", file)
	d = rbind(d, tmp)
}

save(d, file="data.result.age.sites.RData")




###
ggplot(p) +
	#geom_density_2d(aes(time, PostMode * Ne), bins = 5, size = 1/3, color = "black", alpha = 0.8) +
	#geom_linerange(aes(x = x, ymin = prev, ymax = time)) +
	geom_raster(aes(x, PostMode * Ne), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_raster(aes(mid, parent.time), position = "identity", stat = "bin2d", binwidth = c(5/100, 5/100)) +
	#geom_smooth(data = pp, aes(x, Lower * Ne), method = "lm", color = "darkred") +
	#geom_smooth(data = pp, aes(x, Upper * Ne), method = "lm", color = "darkgreen") +
	geom_smooth(data = pp, aes(x, node.time), method = "lm", size = 1/2, color = "grey40") +
	geom_smooth(data = pp, aes(x, mid), method = "lm", size = 1/2, color = "grey40") +
	geom_smooth(data = pp, aes(x, parent.time), method = "lm", size = 1/2, color = "grey40") +
	geom_abline(intercept = c(0,0), slope = 1) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "black", size = 2/3, se = F) +
	geom_smooth(data = pp, aes(x, PostMode * Ne), method = "lm", color = "white", linetype = "22", size = 2/3, fill = "black") +
	geom_segment(data = tm, aes(x=x, xend=xend, y=y, yend=yend), color="grey50") +
	# geom_point(aes(x, Lower * Ne), alpha=0.25, color = "darkred") +
	# geom_point(aes(x, Upper * Ne), alpha=0.25, color = "darkgreen") +
	# geom_point(aes(x, PostMode * Ne), alpha=0.25) +
	#geom_point(aes(time, prev)) +
	coord_cartesian(xlim = c(0.5, 40000), ylim = c(0.5, 40000), expand = F) +
	scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
	#scale_fill_gradientn(colours = c("white", "deepskyblue", "royalblue2", "magenta4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	#scale_fill_gradientn(colours = c("yellow", "orange", "orangered", "darkred", "coral4"), na.value = "navy", trans = "log10", limits = c(1, 1000)) +
	scale_fill_gradientn(colours = c("lightcyan2", "deepskyblue", "blue2", "navy"), na.value = "black", trans = "log10", limits = c(1, 1000)) +
	#scale_color_gradientn(colours = c("lightblue", "deepskyblue", "royalblue", "darkblue"), na.value = "navy", limits = c(0, 4)) +
	#scale_fill_gradient(low = "green", high = "red") +
	scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
	theme_few() +
	theme(aspect.ratio = 1,
				#panel.spacing.y = unit(0.05, "cm"),
				#plot.margin = unit(c(-0.1, 0.05, 0.05, 0.05), "cm"),
				axis.text = element_text(size = 8),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				#panel.background = element_rect(fill = "grey90"),
				strip.text.x = element_text(face = "bold", hjust = 0, size = 8),
				strip.text.y = element_text(face = "bold.italic", angle = 0),
				legend.justification=c(1,0), legend.position=c(0.995,1-0.995), legend.box.background = element_rect(colour = "black"),
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Inferred age") + xlab("True age")
###




###
stop()
###

library(data.table)

files = dir(pattern=".+\\.age\\.sites\\.HMM\\.([a-zA-Z]+)\\.txt$")
files = c(files, dir(pattern=".+\\.age\\.sites\\.HMM\\.([a-zA-Z]+)\\.HardBreaks\\.txt$"))

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	if (grepl("HardBreaks", file)) {
		tmp$type = "Hard breaks"
		tmp$method = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	} else {
		tmp$type = "Soft breaks"
		tmp$method = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.sites\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	}
	d = rbind(d, tmp)
}

save(d, file="data.result.age.sites.HMM.check.RData")



files = dir(pattern=".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$")
files = c(files, dir(pattern=".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$"))

d = NULL

for (file in files) {
	print(file)
	tmp = fread(file, header = T)
	if (grepl("HardBreaks", file)) {
		tmp$type = "Hard breaks"
		tmp$method = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.HardBreaks\\.txt$", "\\2", file)
	} else {
		tmp$type = "Soft breaks"
		tmp$method = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\1", file)
		tmp$clock  = sub(".+\\.age\\.pairs\\.([A-Z]+)\\.([a-zA-Z]+)\\.txt$", "\\2", file)
	}
	d = rbind(d, tmp)
}

save(d, file="data.result.age.pairs.check.RData")







###############


sim = NULL

d = fread("./truth/simres.age.sites.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Mutation clock"
d$method = "True IBD"
sim = rbind(sim, d)

d = fread("./truth/simres.age.sites.SIM.Rec.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Recombination clock"
d$method = "True IBD"
sim = rbind(sim, d)

d = fread("./truth/simres.age.sites.SIM.MutRec.HardBreaks.txt", header = T, stringsAsFactors = F)
d$tag = ""
d$error = NA
d$clock = "Combined clock"
d$method = "True IBD"
sim = rbind(sim, d)


sim = cbind(sim, times[sim$MarkerID + 1, ])


save(sim, file="result.truth.RData")




#############



d = fread("_tmp.txt", header = T, stringsAsFactors = F)

d$a = sprintf("%d.%s", d$SampleID0, d$Chr0)
d$b = sprintf("%d.%s", d$SampleID1, d$Chr1)

d = split(d, d$Shared)

s = d$`1`
o = d$`0`

hh = unique(s$a, s$b)

for (h in hh) {
	
	i = which(o$a == h | o$b == h)
	
	lhs = o$SegmentLHS[i]
	rhs = o$SegmentRHS[i]
	
	j = which(s$a == h | s$b == h)
	
	#if (any(lhs < min(s$SegmentLHS[j]))) cat("L")
	#if (any(rhs > max(s$SegmentRHS[j]))) cat("R")
	
	if (rhs - lhs > max(s$SegmentRHS[j] - s$SegmentLHS[j])) cat(".")
}



mut = as.data.table(read.table("../OutOfAfricaHapMap20.mutations.txt", header = T, stringsAsFactors = F))
rec = as.data.table(read.table("../OutOfAfricaHapMap20.records.txt", header = T, stringsAsFactors = F))
tim = as.data.table(read.table("./OutOfAfricaHapMap20.times.txt", header = T, stringsAsFactors = F))

rec$a = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\1", rec$children))
rec$b = as.numeric(sub("^([0-9]+),([0-9]+)$", "\\2", rec$children))





nn = fread("nn_tru_H_id_100.age.sites.HMM.Rec.HardBreaks.txt", header = T, stringsAsFactors = F)
rd = fread("rd_tru_H_id_100.age.sites.HMM.Rec.HardBreaks.txt", header = T, stringsAsFactors = F)

nn$ci = nn$PostCI975 - nn$PostCI025
rd$ci = rd$PostCI975 - rd$PostCI025

nn$ci = nn$Upper - nn$Lower
rd$ci = rd$Upper - rd$Lower

plot(ecdf(nn$ci * 20000)(seq(0, 20000, length.out = 500)), type = 'l', col = "blue")
lines(ecdf(rd$ci * 20000)(seq(0, 20000, length.out = 500)), type = 'l', col="red")




sim = fread("~/Research/vanilla/simres.age.pairs.SIM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
sim = fread("dev_err_H_id_100_NN.age.pairs.HMM.Mut.HardBreaks.txt", header = T, stringsAsFactors = F)
sim = fread("~/Research/rvage/africa/test.age.pairs.HMM.MutRec.HardBreaks.txt", header = T, stringsAsFactors = F)

sim = split(sim, sim$MarkerID)

res = c()

for (tag in names(sim)) {
	x = sim[[tag]]
	s = which(x$Shared == 1)
	o = which(x$Shared == 0)
	shr = x[s,]
	oth = x[o,]
	
	if (any(shr$SampleID0 == shr$SampleID1)) next
	if (any(oth$SampleID0 == oth$SampleID1)) next
	
	sam = unique(c(shr$SampleID0, shr$SampleID1))
	
	lhs = list()
	rhs = list()
	for (s in sam) {
		i = which(shr$SampleID0 == s | shr$SampleID1 == s)
		lhs[[ as.character(s) ]] = min(shr$SegmentLHS[i])
		rhs[[ as.character(s) ]] = max(shr$SegmentRHS[i])
	}
	
	unr = unique(c(oth$SampleID0, oth$SampleID1))
	unr = setdiff(unr, sam)
	
	flag = F
	for (u in unr) {
		i = which(oth$SampleID0 == u | oth$SampleID1 == u)
		s = intersect(sam, c(oth$SampleID0[i], oth$SampleID1[i]))
		wl = which(oth$SegmentLHS[i] < min(unlist(lhs[as.character(s)])))
		wr = which(oth$SegmentRHS[i] > max(unlist(rhs[as.character(s)])))
		n = length(unique(c(wl, wr)))
		if (n > 0) {
			cat(n, "")
			flag = T
		}
	}
	if (flag)
		cat("  of", length(unr), " at ", length(sam), "\n")
}
range(res)



