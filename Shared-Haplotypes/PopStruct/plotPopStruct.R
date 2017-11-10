
library(ggplot2)
library(reshape2)

colour.grp = c(AFR="#FFCD33",
							 "AFR/AMR"="#FF9900",
							 AMR="#FF3D3D",
							 EAS="#ADFF33",
							 EUR="#64EBFF",
							 SAS="#FF30FF")

colour.pop = c(ESN="#FFCD00",
							 GWD="#FFB900",
							 LWK="#CC9933",
							 MSL="#E1B919",
							 YRI="#FFB933",
							 ACB="#FF9900",
							 ASW="#FF6600",
							 CLM="#CC3333",
							 MXL="#E10033",
							 PEL="#FF0000",
							 PUR="#CC3300",
							 CDX="#339900",
							 CHB="#ADCD00",
							 CHS="#00FF00",
							 JPT="#008B00",
							 KHV="#00CC33",
							 CEU="#0000FF",
							 FIN="#00C5CD",
							 GBR="#00EBFF",
							 IBS="#6495ED",
							 TSI="#00008B",
							 BEB="#8B008B",
							 GIH="#9400D3",
							 ITU="#B03060",
							 PJL="#E11289",
							 STU="#FF00FF")



panel = read.table("integrated_call_samples_v3.20130502.ALL.panel", header = TRUE, stringsAsFactors = FALSE)
panel$super_pop[which(panel$pop == "ACB")] = "AFR/AMR"
panel$super_pop[which(panel$pop == "ASW")] = "AFR/AMR"
panel = panel[order(panel$sample), ]
panel = panel[order(panel$pop), ]
panel = panel[order(panel$super_pop), ]
rownames(panel) = panel$sample




files = dir(pattern = "^chr([0-9]+)\\.f05\\.shared$")
m = list()
for (file in files) {
	m[[file]] = as.matrix(read.table(file, header = TRUE, row.names = ".", stringsAsFactors = FALSE))
}
M = matrix(0, nrow(m[[1]]), ncol(m[[1]]), dimnames = dimnames(m[[1]]))
for (i in 1:length(m)) {
	M = M + m[[i]]
}
F05 = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
for (col in panel$sample) {
	F05[, col] = M[panel$sample, col]
}



files = dir(pattern = "^chr([0-9]+)\\.f25\\.shared$")
m = list()
for (file in files) {
	m[[file]] = as.matrix(read.table(file, header = TRUE, row.names = ".", stringsAsFactors = FALSE))
}
M = matrix(0, nrow(m[[1]]), ncol(m[[1]]), dimnames = dimnames(m[[1]]))
for (i in 1:length(m)) {
	M = M + m[[i]]
}
F25 = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
for (col in panel$sample) {
	F25[, col] = M[panel$sample, col]
}



files = dir(pattern = "^chr([0-9]+)\\.f50\\.shared$")
m = list()
for (file in files) {
	m[[file]] = as.matrix(read.table(file, header = TRUE, row.names = ".", stringsAsFactors = FALSE))
}
M = matrix(0, nrow(m[[1]]), ncol(m[[1]]), dimnames = dimnames(m[[1]]))
for (i in 1:length(m)) {
	M = M + m[[i]]
}
F50 = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
for (col in panel$sample) {
	F50[, col] = M[panel$sample, col]
}



files = dir(pattern = "^chr([0-9]+)\\.f100\\.shared$")
m = list()
for (file in files) {
	m[[file]] = as.matrix(read.table(file, header = TRUE, row.names = ".", stringsAsFactors = FALSE))
}
M = matrix(0, nrow(m[[1]]), ncol(m[[1]]), dimnames = dimnames(m[[1]]))
for (i in 1:length(m)) {
	M = M + m[[i]]
}
F100 = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
for (col in panel$sample) {
	F100[, col] = M[panel$sample, col]
}


F00 = F100-F05

X = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
X[upper.tri(F05)] = F05[upper.tri(F05)]
X[lower.tri(F25)] = F25[lower.tri(F25)]


X = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
X[upper.tri(F50)] = F50[upper.tri(F50)]
X[lower.tri(F100)] = F100[lower.tri(F100)]


X = matrix(0, nrow(M), ncol(M), dimnames = list(panel$sample, panel$sample))
X[upper.tri(F05)] = F05[upper.tri(F05)]
X[lower.tri(F00)] = F00[lower.tri(F00)]



d = melt(X, varnames = c("a", "b"))
d = cbind(d, a.pop = panel$pop[d$a], a.grp = panel$super_pop[d$a])
d = cbind(d, b.pop = panel$pop[d$b], b.grp = panel$super_pop[d$b])


tmp = as.numeric(as.factor(d$a))

d.pop = NULL
for (pop in unique(panel$pop)) {
	m = mean(range(tmp[which(d$a.pop == pop)]))
	d.pop = rbind(d.pop, data.frame(pop = pop, m = m))
}

d.grp = NULL
for (grp in unique(panel$super_pop)) {
	m = mean(range(tmp[which(d$a.grp == grp)]))
	d.grp = rbind(d.grp, data.frame(grp = grp, m = m))
}


gg = ggplot(data=d) + 
	geom_raster(aes(x=a, y=b, fill=log(value+1))) + 
	geom_segment(aes(x=a, xend=a, y=-5,  yend=-55,  colour=a.pop)) +
	geom_segment(aes(x=a, xend=a, y=-60, yend=-110, colour=a.grp)) +
	geom_segment(aes(y=b, yend=b, x=-5,  xend=-55,  colour=b.pop)) +
	geom_segment(aes(y=b, yend=b, x=-60, xend=-110, colour=b.grp))

for(theta in seq(pi/8, 2*pi, length.out=16)) {
	gg = gg + 
		geom_text(data = d.pop, bquote(aes(x=m+.(cos(theta)*2), y=-30+.(sin(theta)*2), label=pop)), colour="black", alpha=1/8) +
		geom_text(data = d.pop, bquote(aes(y=m+.(cos(theta)*2), x=-30+.(sin(theta)*2), label=pop, angle=-90)), colour="black", alpha=1/8) +
		geom_text(data = d.grp, bquote(aes(x=m+.(cos(theta)*2), y=-85+.(sin(theta)*2), label=grp)), fontface="bold", colour="white", alpha=1/8) +
		geom_text(data = d.grp, bquote(aes(y=m+.(cos(theta)*2), x=-85+.(sin(theta)*2), label=grp, angle=-90)), fontface="bold", colour="white", alpha=1/8)
}

gg = gg +
	geom_text(data = d.pop, aes(x=m, y=-30, label=pop), colour="white") +
	geom_text(data = d.pop, aes(y=m, x=-30, label=pop, angle=-90), colour="white") +
	geom_text(data = d.grp, aes(x=m, y=-85, label=grp), fontface="bold", colour="black") +
	geom_text(data = d.grp, aes(y=m, x=-85, label=grp, angle=-90), fontface="bold", colour="black") +
	scale_fill_gradient(low = "white", high = "black") +
	scale_colour_manual(values=c(colour.pop, colour.grp)) +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	coord_fixed() +
	theme(axis.title=element_blank(),
				axis.text=element_blank(),
				axis.ticks=element_blank(),
				axis.line=element_blank()) #),
				#legend.position="none")

ggsave(gg, filename = "plot.popstruct.png", width = 17, height = 15)





