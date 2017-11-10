

library(data.table)
library(ggplot2)
library(ggthemes)

se <- function(x) sqrt(var(x)/length(x))



load("../data.OutOfAfricaHapMap20.G.RData")
A=G
load("./data.OutOfAfricaHapMap20.GenErr_1000G.G.RData")
B=G

rm(G)

fa = rowSums(A)
fb = rowSums(B)

ia = which(fa >= 1 & fa <= 250)
ib = which(fb >= 1 & fb <= 250)

ga = split(ia, fa[ia])
gb = split(ib, fb[ib])


fn = list()

for (k in names(ga)) {
	print(k)
	
	aa = A[ga[[k]], ]
	bb = B[ga[[k]], ]
	
	x = which(aa == 1)
	
	fn[[k]] = table(bb[x])
}


fp = list()

for (k in names(gb)) {
	print(k)
	
	aa = A[gb[[k]], ]
	bb = B[gb[[k]], ]
	
	x = which(bb == 1)
	
	fp[[k]] = table(aa[x])
}


d = NULL

for (k in names(ga)) {
	
	n = fn[[k]]
	
	n0 = 0
	n2 = 0
	
	if ("0" %in% names(n)) n0 = n["0"]
	if ("2" %in% names(n)) n2 = n["2"]
	
	p = fp[[k]]
	
	p0 = 0
	p2 = 0
	
	if ("0" %in% names(p)) p0 = p["0"]
	if ("2" %in% names(p)) p2 = p["2"]
	
	d = rbind(d, 
						data.table(fk = k, type = "False negative", mode = "0", n = n0 / sum(n)),
						data.table(fk = k, type = "False negative", mode = "2", n = n2 / sum(n)),
						data.table(fk = k, type = "False positive", mode = "0", n = p0 / sum(p)),
						data.table(fk = k, type = "False positive", mode = "2", n = p2 / sum(p)))
		
}

d$fk = factor(d$fk, levels = names(ga), ordered = T)


fpn = ggplot() +
	geom_bar(data = d[which(d$type == "False positive"), ], aes(x=fk, y=n, fill = mode), stat = "identity", color = "grey60", show.legend = F) +
	geom_bar(data = d[which(d$type == "False negative"), ], aes(x=fk, y=n*-1, fill = mode), stat = "identity", color="grey60", show.legend = F) +
	geom_hline(yintercept = 0, size = 0.5) +
	scale_x_discrete(breaks = c(1, seq(10, 500, by=10)), labels = c(1, seq(10, 500, by=10))/50) +
	scale_y_continuous(breaks = seq(-20, 20, by = 2)/100, labels = abs(seq(-20, 20, by = 2))/100) +
	scale_fill_manual(values = c("grey60", "grey60")) +
	theme_few() +
	theme(#aspect.ratio=1/2,
				panel.border=element_rect(fill = NA, colour = "black", size = 2/3),
				legend.title=element_blank(),
				panel.grid = element_line(),
				panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank()) +
	ylab("    FNR                     FPR") +
	xlab("Allele frequency (%)")
fpn

ggsave(fpn, filename = "___plot.falseposneg.pdf", width = 9, height = 3)


mean(d$n[which(d$type == "False positive")]);  se(d$n[which(d$type == "False positive")]);  
mean(d$n[which(d$type == "False negative")]);  se(d$n[which(d$type == "False negative")]);  


p = NULL

for (k in 1:250) {
	x = which(fa == k)
	y = table(fb[x])
	y = y / sum(y)
	
	p = rbind(p, data.table(fk = k, mk = as.numeric(names(y)), p = as.vector(y)))
}

# del = which(p$fk == p$mk)
# p = p[-del,]

fkdiff = ggplot(p) +
	geom_raster(aes(fk, mk-fk, fill = p)) +
	scale_fill_gradientn(colours = c("grey", "tomato2", "darkred", "red4", "black"), na.value = "black", breaks = (1:9)/10, limits = c(0, 1)) +
	scale_x_continuous(breaks = c(1, seq(10, 500, by=10)), labels = c(1, seq(10, 500, by=10))/50, minor_breaks = (1:250)-0.5) +
	scale_y_continuous(minor_breaks = (-50:50)-0.5) +
	coord_cartesian(xlim = c(0.01, 5.01)*50, ylim=c(-10.5, 10.5), expand = F) +
	theme_few() +
	theme(#aspect.ratio=1/2,
		panel.border=element_rect(fill = NA, colour = "black", size = 2/3),
		legend.title=element_blank(),
		legend.key.height = unit(0.25, "cm"),
		legend.key.width = unit(1.25, "cm"),
		legend.margin = margin(0.1,0.1,0.1,0, unit = "cm"),
		legend.background = element_rect(fill = "white", colour = "black", size = 1/3),
		legend.direction = "horizontal",
		legend.position = c(0.995, 0.015), legend.justification = c(1, 0),
		legend.text = element_text(size = 8),
		panel.grid = element_line(),
		panel.grid.major = element_blank(),
		panel.grid.minor.y = element_line(size = 0.5, colour = "grey90"),
		panel.grid.minor.x = element_line(size = 0.5, colour = "grey90"),
		panel.grid.minor = element_line(size = 0.5, colour = "grey90")) +
	ylab("Allele count difference") +
	xlab("Allele frequency (%)")
fkdiff

ggsave(fkdiff, filename = "___plot.fk_diff.pdf", width = 9, height = 3)



