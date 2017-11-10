
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)


# files = dir(pattern = "^count\\.chr[0-9]+\\.psm\\.[0-9]+\\.txt$")
# 
# d = fread(files[1], header = T, stringsAsFactors = F)
# d$N = 0
# 
# for (file in files) {
# 	cat(".")
# 	
# 	tmp = fread(file, header = T, stringsAsFactors = F)
# 	
# 	if (nrow(d) != nrow(tmp)) {
# 		print(file)
# 		next
# 	}
# 	
# 	d$N = d$N + tmp$N
# }
# 
# n = d
# rm(d)
# save(n, file = "_data.count.psm.RData")
# 
# #####
# 
# load("_data.count.psm.RData")


#####

load("data.psm.RData")









panel = read.table("../data/integrated_call_samples_v3.20130502.ALL.panel", header = TRUE, stringsAsFactors = FALSE)
panel$super_pop[which(panel$pop == "ACB")] = "AFR/AMR"
panel$super_pop[which(panel$pop == "ASW")] = "AFR/AMR"
panel = panel[order(panel$sample), ]
panel = panel[order(panel$pop), ]
panel = panel[order(panel$super_pop), ]
rownames(panel) = panel$sample


pop = panel$pop
names(pop) = panel$sample

sup = panel$super_pop
names(sup) = panel$sample

ids = 1:nrow(panel)
names(ids) = panel$sample



m = matrix(NA, nrow(panel), nrow(panel), dimnames = list(panel$sample, panel$sample))



sample = read.table("../data/1000G.chr1_c.sample.txt", header = TRUE, stringsAsFactors = FALSE)

n$id0 = sample$Label[ n$SampleID0 + 1 ]
n$id1 = sample$Label[ n$SampleID1 + 1 ]


coord = as.matrix(data.frame(x = n$id0, y = n$id1))
m[ coord ] = n$N

coord = as.matrix(data.frame(x = n$id1, y = n$id0))
m[ coord ] = n$N

if (any(m == 0))
	m[which(m == 0)] = NA


save(m, file = "_data.matrix.psm.RData")













n$pop0 = pop[ n$id0 ]
n$pop1 = pop[ n$id1 ]
n$sup0 = sup[ n$id0 ]
n$sup1 = sup[ n$id1 ]
n$ord0 = ids[ n$id0 ]
n$ord1 = ids[ n$id1 ]

n$id0 = as.factor(n$id0)
n$id1 = as.factor(n$id1)
n$pop0 = as.factor(n$pop0)
n$pop1 = as.factor(n$pop1)
n$sup0 = as.factor(n$sup0)
n$sup1 = as.factor(n$sup1)
n$ord0 = as.factor(n$ord0)
n$ord1 = as.factor(n$ord1)


share = n
rm(n)

save(share, file = "_data.share.psm.RData")





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




share$P0 = as.factor(sprintf("%s %s", as.character(share$sup0), as.character(share$pop0)))
share$P1 = as.factor(sprintf("%s %s", as.character(share$sup1), as.character(share$pop1)))


z = sapply(split(share, share$sup0), function(x) sapply(split(x, x$sup1), function(y) sum(y$N)))
if (any(is.na(z))) z[which(is.na(z))] = 0
z = t(z)
z = z/rowSums(z)*100


z = sapply(split(share, share$P0), function(x) sapply(split(x, x$P1), function(y) sum(y$N)/nrow(y)))
if (any(is.na(z))) z[which(is.na(z))] = 0
z = t(z)
z = z/rowSums(z)*100

m = array(0, dim = dim(z), dimnames = dimnames(z))
for (i in colnames(z)) {
	for (j in colnames(z)) {
		if (i == j) 	
			m[i,j] = z[i, j]
		else
			m[i,j] = z[i, j] + z[j, i]
	}
}
m[lower.tri(m)] = NA

ggplot(melt(z)) + geom_raster(aes(Var1, Var2, fill = value))


gg = ggplot(data=share) + 
	geom_raster(aes(x=ord0, y=ord1, fill=log10(N+1))) + 
	geom_segment(aes(x=ord0, xend=ord0, y=-5,  yend=-55,  colour=pop0)) +
	geom_segment(aes(x=ord0, xend=ord0, y=-60, yend=-110, colour=sup0)) +
	geom_segment(aes(y=ord1, yend=ord1, x=2504+5,  xend=2504+55,  colour=pop1)) +
	geom_segment(aes(y=ord1, yend=ord1, x=2504+60, xend=2504+110, colour=sup1)) +
	# geom_text(data = d.pop, aes(x=m, y=-30, label=pop), colour="white") +
	# geom_text(data = d.pop, aes(y=m, x=-30, label=pop, angle=-90), colour="white") +
	# geom_text(data = d.grp, aes(x=m, y=-85, label=grp), fontface="bold", colour="black") +
	# geom_text(data = d.grp, aes(y=m, x=-85, label=grp, angle=-90), fontface="bold", colour="black") +
	scale_fill_gradient(low = "white", high = "black") +
	scale_colour_manual(values=c(colour.pop, colour.grp)) +
	scale_x_discrete(expand = c(0,0)) +
	scale_y_discrete(expand = c(0,0)) +
	theme_classic() +
	theme(aspect.ratio = 1,
				axis.title=element_blank(),
				axis.text=element_blank(),
				axis.ticks=element_blank(),
				axis.line=element_blank())






