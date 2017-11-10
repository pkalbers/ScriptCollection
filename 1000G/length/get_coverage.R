

library(data.table)
library(ggplot2)
library(ggthemes)


args = commandArgs(T)

ibd = fread(args[1], header = T, stringsAsFactors = F)


### unique segments per target individual


n = max(ibd$RHS)


### DGT

stat = NULL # collect stats
segm = list()

tars = sort(unique(c(ibd$SampleID0, ibd$SampleID1)))

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
	
	brk = which(lhs != 0)
	lhs[brk] = lhs[brk] + 1
	
	brk = which(rhs != n)
	rhs[brk] = rhs[brk] - 1
	
	d = data.frame(fk = ibd$Fk[x], lhs = lhs, rhs = rhs, key = sprintf("%d %d", lhs, rhs))
	d = d[order(d$fk, d$key), ]
	
	del = which(duplicated(d$key))
	if (length(del) > 0) d = d[-del, ]
	
	d$key = NULL
	segm[[ as.character(tar) ]] = d
}



cover = list()

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
		
		if (sum(sapply(rng, length)) == length(unlist(rng))) {
			
			rng = unlist(rng)
			c[rng] = c[rng] + 1
			
		} else {
			
			for (i in 1:length(rng)) {
				sub = rng[[i]]
				c[sub] = c[sub] + 1
			}
			
		}
		
		tt = data.table(target = as.numeric(tar), fk = as.numeric(k))
		for (f in 1:24) {
			x = length(which(c >= f))
			tt[[ sprintf("fold%d", f) ]] = x / n
		}
		
		tmp = rbind(tmp, tt)
	}
	cat("\n")
	
	cover[[tar]] = tmp
}



save(stat, cover, file = sprintf("cover.%s.RData", agrs[1]))



stop("DONE")





#########

load("stat_dgt.coverage.RData")


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

d = NULL


load("stat_dgt.coverage.RData")

del = which(sapply(cover, nrow) != nrow(cover[[1]]))
if (length(del) > 0) cover[del] = NULL

for (f in c(1:24)) {
	fold = lapply(cover, function(x, f)  as.data.table(t(data.table(x[[ sprintf("fold%d", f) ]]))), f)
	fold = rbindlist(fold)
	d = rbind(d, data.table(fk = 2:25, mn = apply(fold, 2, mean), se = apply(fold, 2, se), fold = sprintf("% 3d-fold", f), type = "(b) Detected"))
}


load("stat_tru.coverage.RData")

del = which(sapply(cover, nrow) != nrow(cover[[1]]))
if (length(del) > 0) cover[del] = NULL

for (f in c(1:24)) {
	fold = lapply(cover, function(x, f)  as.data.table(t(data.table(x[[ sprintf("fold%d", f) ]]))), f)
	fold = rbindlist(fold)
	d = rbind(d, data.table(fk = 2:25, mn = apply(fold, 2, mean), se = apply(fold, 2, se), fold = sprintf("% 3d-fold", f), type = "(a) Truth"))
}





bin = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = -0.1, y1 = 1.1)
del = which(1:nrow(bin) %% 2 == 0)
bin = bin[-del, ]

cols = colorRampPalette(c(colorRampPalette(c("grey70", "pink"))(3)[2], 
													colorRampPalette(c("grey70", "blue"))(3)[2], 
													colorRampPalette(c("grey20", "red"))(3)[2], 
													colorRampPalette(c("grey90", "darkorange"))(3)[2]))(24)
cols = rev(cols)

gg = ggplot(d) +
	facet_grid(.~type) +
	geom_rect(data = bin, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1), fill = "grey90", alpha = 0.5) +
	#geom_ribbon(aes(x = fk, ymin = mn-se, ymax = mn+se, fill = fold), alpha = 0.5) +
	geom_line(aes(x = fk, y = mn, color = fold)) +
	geom_point(aes(x = fk, y = mn, color = fold, size = se), alpha = 0.5) +
	geom_point(aes(x = fk, y = mn, color = fold)) +
	scale_x_continuous(breaks = 2:25) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
	scale_color_manual(values = cols) +
	coord_cartesian(xlim = c(1.5, 25.5), ylim = c(-0.01, 1.01), expand = F) +
	theme_few() +
	theme(aspect.ratio = 9/8, 
				strip.text = element_text(face = "bold", hjust = 0),
				legend.key.size = unit(0.6, "cm"),
				legend.text = element_text(size = 8),
				legend.title = element_text(size = 9),
				legend.key.height = unit(0.4, "cm"),
				panel.border = element_rect(fill = NA, colour = "black", size = 1/2),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") + ylab("Mean cumulative coverage (~SE)") + 
	guides(size = guide_legend(title = "SE"), color = guide_legend(title = "Coverage", ncol = 1, override.aes = list(size = 1)))
gg

ggsave(gg, filename = "_plot.phase_coverage.pdf", width = 11, height = 6)



