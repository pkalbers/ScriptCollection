

library(data.table)
library(ggplot2)
library(ggthemes)



mrk = fread("1000G.chr20.marker.txt", header = T, stringsAsFactors = F)
ibd = fread("result.chr20.ibd.FGT.txt", header = T, stringsAsFactors = F)
ibd = fread("result.chr20.ibd.DGT.txt", header = T, stringsAsFactors = F)
ibd = fread("result.chr20.ibd.HMM.txt", header = T, stringsAsFactors = F)


POS = mrk$Position


folds = c(1:5, seq(10, 50, by =10))

tars = sort(unique(c(ibd$SampleID0, ibd$SampleID1)))
tars = sample(tars, 500)


### unique segments per target individual



### ANY!!!

stat = NULL # collect stats
segm = list()

for (tar in tars) {
	print(tar)
	
	x = unique(c(which(ibd$SampleID0 == tar), which(ibd$SampleID1 == tar)))
	
	n.c = length(unique(c(ibd$SampleID0[x], ibd$SampleID1[x]))) - 1
	
	d = data.frame(lhs = ibd$LHS[x], rhs = ibd$RHS[x])
	
	n.a = nrow(d)
	
	d = unique(d)
	
	n.b = nrow(d)
	
	stat = rbind(stat, data.table(target = tar, nothers = n.c, npair = n.a, usegm = n.b))
	
	
	lhs = ibd$LHS[x]
	rhs = ibd$RHS[x]
	
	oth = data.table(a=ibd$SampleID0[x], b=ibd$SampleID1[x])
	r = which(oth$a == tar); if (length(r) > 0) oth$a[r] = 0
	r = which(oth$b == tar); if (length(r) > 0) oth$b[r] = 0
	oth = rowSums(oth)
	
	d = data.frame(fk = ibd$Fk[x], lhs = lhs, rhs = rhs, key = sprintf("%d %d %d", lhs, rhs, oth))
	d = d[order(d$fk, d$key), ]
	
	del = which(duplicated(d$key))
	if (length(del) > 0) d = d[-del, ]
	
	d$key = NULL
	segm[[ as.character(tar) ]] = d
}

cover = list()

n = length(POS)

for (tar in names(segm)) {
	print(tar)
	
	c = rep(0, n)
	d = segm[[tar]]
	d = split(d, d$fk)
	
	tmp = NULL
	
	for (k in names(d)) {
		cat(".")
		rng = matrix(c(d[[k]]$lhs, d[[k]]$rhs), ncol = 2, byrow = F)
		rng = apply(rng, 1, function(x) (x[1]):(x[2]))
		
		for (i in 1:length(rng)) {
			sub = rng[[i]]
			c[sub] = c[sub] + 1
		}
		
		tt = data.table(target = as.numeric(tar), fk = as.numeric(k))
		for (f in folds) {
			x = length(which(c >= f))
			tt[[ sprintf("fold%d", f) ]] = x / n
		}
		
		tmp = rbind(tmp, tt)
	}
	cat("\n")
	
	cover[[tar]] = tmp
}

save(stat, cover, file = "stat_dgt.coverage.RData")
save(stat, cover, file = "stat_fgt.coverage.RData")
save(stat, cover, file = "stat_hmm.coverage.RData")






stop()


#########

load("stat_dgt.coverage.RData")
load("../../1000G/ibd/stat_dgt.coverage.RData")


se <- function(x) sqrt(var(x)/length(x))


# mean N others
mean(stat$nothers)
se(stat$nothers)

plot(density(stat$nothers))

# mean N segments
mean(stat$npair)
se(stat$npair)

plot(density(stat$npair))

# mean N unique segments
mean(stat$usegm)
se(stat$usegm)

plot(density(stat$usegm))


# percent FOLD1 not reached
length(which(!sapply(cover, function(x) (any(x$fold1 == 1)) ))) / length(cover) * 100



### mean coverage

se <- function(x) sqrt(var(x)/length(x))

d = NULL


load("stat_fgt.coverage.RData")

del = which(sapply(cover, nrow) != nrow(cover[[1]]))
if (length(del) > 0) cover[del] = NULL

for (f in folds) {
	fold = lapply(cover, function(x, f)  as.data.table(t(data.table(x[[ sprintf("fold%d", f) ]]))), f)
	fold = rbindlist(fold)
	d = rbind(d, data.table(fk = 2:25, mn = apply(fold, 2, mean), se = apply(fold, 2, se), fold = sprintf("% 3d-fold", f), type = "(a) FGT, phased haplotypes"))
}

load("stat_dgt.coverage.RData")

del = which(sapply(cover, nrow) != nrow(cover[[1]]))
if (length(del) > 0) cover[del] = NULL

for (f in folds) {
	fold = lapply(cover, function(x, f)  as.data.table(t(data.table(x[[ sprintf("fold%d", f) ]]))), f)
	fold = rbindlist(fold)
	d = rbind(d, data.table(fk = 2:25, mn = apply(fold, 2, mean), se = apply(fold, 2, se), fold = sprintf("% 3d-fold", f), type = "(b) DGT, genotypes"))
}




load("stat_hmm.coverage.RData")

del = which(sapply(cover, nrow) != nrow(cover[[1]]))
if (length(del) > 0) cover[del] = NULL

for (f in folds) {
	fold = lapply(cover, function(x, f)  as.data.table(t(data.table(x[[ sprintf("fold%d", f) ]]))), f)
	fold = rbindlist(fold)
	d = rbind(d, data.table(fk = 2:25, mn = apply(fold, 2, mean), se = apply(fold, 2, se), fold = sprintf("% 3d-fold", f), type = "HMM, genotypes"))
}



bin = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = -0.1, y1 = 1.1)
del = which(1:nrow(bin) %% 2 == 0)
bin = bin[-del, ]

cols = colorRampPalette(c(colorRampPalette(c("grey70", "pink2"))(3)[2], 
													colorRampPalette(c("grey70", "blue"))(3)[2], 
													colorRampPalette(c("grey20", "red"))(3)[2], 
													colorRampPalette(c("grey90", "darkorange"))(3)[2]))(10)
cols = rev(cols)

gg = ggplot(d) +
	facet_wrap(~type) +
	geom_rect(data = bin, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1), fill = "grey90", alpha = 0.5) +
	#geom_ribbon(aes(x = fk, ymin = mn-se, ymax = mn+se, fill = fold), alpha = 0.5) +
	geom_line(aes(x = fk, y = mn, color = fold)) +
	geom_point(aes(x = fk, y = mn, color = fold, size = se), alpha = 0.5) +
	geom_point(aes(x = fk, y = mn, color = fold), size = 0.5, show.legend = F) +
	scale_x_continuous(breaks = 2:25, labels = c('2','','','5','','','','','10','','','','','15','','','','','20','','','','','25')) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
	scale_color_manual(values = cols) +
	coord_cartesian(xlim = c(1.5, 25.5), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 1, 
				strip.text = element_text(face = "bold", hjust = 0),
				#axis.text.x = element_text(size = 9),
				legend.key.size = unit(0.6, "cm"),
				legend.text = element_text(size = 9),
				legend.title = element_text(size = 9),
				#legend.key.height = unit(0.6, "cm"),
				panel.border = element_rect(fill = NA, colour = "black", size = 1/2),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") + ylab("Mean chromosome coverage (~SE)") + 
	guides(size = guide_legend(title = "SE"), color = guide_legend(title = "Coverage", ncol = 1, label.hjust = 1, override.aes = list(size = 1.5)))
gg

ggsave(gg, filename = "_plot.phase_coverage_1kg_hmm.pdf", width = 10, height = 5)









p = match

p$key = sprintf("%d %d %.3f %.3f", p$g0, p$g1, p$g.lhs, p$g.rhs)
del = which(duplicated(p$key))
p = p[-del,]

t = 123

k = which(p$g0 == t | p$g1 == t)
