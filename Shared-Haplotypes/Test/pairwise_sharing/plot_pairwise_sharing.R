#
# plot pairwise sharing matrix
#


library(ggplot2)
library(colorRamps)



split.col <- function(col, as_numeric = TRUE) {
	col = strsplit(col, "|", TRUE)
	
	if (as_numeric) {
		col = lapply(col, as.numeric)
	}
	
	col = Reduce(rbind, col)
	rownames(col) <- NULL
	col
}



plot.pairws <- function(file) {
	
	cat(sprintf("%s\n", file))
	
	target <- sub(".+__(.+)__(.+)__(.+)\\..+", "\\1 \\2 \\3", file);
	target <- paste(rev(strsplit(target, " ", T)[[1]]), collapse = " ")
	
	tab <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
	
	tab <- cbind(tab,
							 sid = as.numeric(sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\1", tab$sample)),
							 key = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\2", tab$sample),
							 pop = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\3", tab$sample),
							 grp = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\4", tab$sample),
							 stringsAsFactors = FALSE)
	
	pmin = min(split.col(tab$lpos))
	pmax = max(split.col(tab$rpos))
	
	xseq <- seq(round(pmin/1000000)*1000000, round(pmax/1000000)*1000000, by = 1000000)
	xseq <- xseq[-1]
	
	tab <- cbind(tab, 
							 index = paste(tab$grp, tab$pop, tab$key), 
							 stringsAsFactors = FALSE)
	
	tab <- tab[order(tab$mpos), ]
	tab <- cbind(tab, order = 1:nrow(tab))
	tab <- tab[order(tab$index), ]
	
	
	plot.blocks <- function(tab) {
		coord.l <- split.col(tab$lpos)
		coord.r <- split.col(tab$rpos)
		
		index <- rev(1:ncol(coord.l))
		alpha <- rev(log(index + 1) / sum(log(index + 1)))
		
		blocks <- function(tab, crd.l, crd.r) {
			rbind(data.frame(name = tab$key, uni = 1:nrow(tab), ymax = tab$order + 0.45, ymin = tab$order - 0.45, xmin = crd.l, xmax = crd.r),
						data.frame(name = tab$key, uni = 1:nrow(tab), ymax = tab$order + 0.45, ymin = tab$order - 0.45, xmin = crd.l, xmax = crd.r))
		}
		
		g = geom_rect(data = blocks(tab, coord.l[, index[1]], coord.r[, index[1]]), aes(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, group = uni), fill = "white")
		
		for (i in 1:length(index)) {
			g = c(g, geom_rect(data = blocks(tab, coord.l[, index[i]], coord.r[, index[i]]), aes(xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, group = uni, fill = name), alpha = alpha[i]))
		}
		
		g
	}
	
	
	plot.breaks <- function(tab) {
		coord.l <- split.col(tab$lpos)
		coord.r <- split.col(tab$rpos)
		
		index <- rev(1:ncol(coord.l))
		
		m <- rbind(data.frame(uni = 1:nrow(tab), x = tab$mpos, y = tab$order + 0.45), 
							 data.frame(uni = 1:nrow(tab), x = tab$mpos, y = tab$order - 0.45))
		
		g = geom_line(data = m, aes(y = y, x = x, group=uni), color="black", size=0.25)
		
		for (i in 1:length(index)) {
			d = rbind(data.frame(uni = 1:nrow(tab), side = paste("l", 1:nrow(tab)), x = coord.l[, index[i]], y = tab$order + 0.45), 
								data.frame(uni = 1:nrow(tab), side = paste("l", 1:nrow(tab)), x = coord.l[, index[i]], y = tab$order - 0.45),
								data.frame(uni = 1:nrow(tab), side = paste("r", 1:nrow(tab)), x = coord.r[, index[i]], y = tab$order + 0.45), 
								data.frame(uni = 1:nrow(tab), side = paste("r", 1:nrow(tab)), x = coord.r[, index[i]], y = tab$order - 0.45))
			g = c(g, geom_line(data = d, aes(y = y, x = x, group=side), color="white", size=0.25))
		}
		
		g = c(g, geom_line(data = m, aes(y = y, x = x, group=uni), color="black", size=0.25))
		
		g
	}
	
	
	hline <- unique(data.frame(y = c(-5, tab$order, max(tab$order) + 5), grp = c(strsplit(target, " ", T)[[1]][1], tab$grp, "")))
	
	
	mline <- rbind(data.frame(uni = 1:nrow(tab), x = tab$mpos, y = tab$order + 0.45), 
								 data.frame(uni = 1:nrow(tab), x = tab$mpos, y = tab$order - 0.45),
								 data.frame(uni = unique(tab$mpos) * -1, x = unique(tab$mpos), y = -5 + 1), 
								 data.frame(uni = unique(tab$mpos) * -1, x = unique(tab$mpos), y = -5 - 1))
	
	lab <- unique(data.frame(name = c(target, tab$index), i = c(-5, tab$order), stringsAsFactors = FALSE))
	
	
	g <- ggplot() + #facet_wrap(~f, ncol = 1) +
		theme_classic() + 
		theme(legend.position="none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 5, hjust = 0)) +
		
		geom_hline(data = hline, aes(yintercept = y), size = 0.5, colour="grey95") +
		
		plot.blocks(tab) + 
		plot.breaks(tab) +
		
		geom_line(data = mline, aes(y = y, x = x, group=uni), color="black", size=0.25) +
		
		scale_x_continuous(expand = c(0, 0), breaks = xseq, labels = xseq / 1000000) +
		scale_y_reverse(expand = c(0, 0), breaks = lab$i, labels = lab$name) +
		xlab("Physical position (Mb)") + ylab("")
	
	outfile <- sprintf("plot_pairws.%s.pdf", file)
	
	ggsave(filename = outfile, width = 50, height = (nrow(hline) + 10) * 0.08, dpi = 300, limitsize = FALSE)
	
}


overlap <- function(tab) {
	coord.l <- split.col(tab$lpos)
	coord.r <- split.col(tab$rpos)
	
	index <- rev(1:ncol(coord.l))
	
	pmin = min(coord.l)
	pmax = max(coord.r)
	
	overl <- NULL
	
	for (i in index) {
		crd.l <- sort(coord.l[, i])
		crd.r <- sort(coord.r[, i])
		
		crd = data.frame(ihg = i, cov = 0, pos = sort(unique(c(coord.l[, i], coord.l[, i] - 1, 
																													 coord.r[, i], coord.r[, i] + 1,
																													 tab$mpos, pmin, pmax))))
		
		blk = data.frame(l = coord.l[, i], m = tab$mpos, r = coord.r[, i])
		blk = blk[order(blk$m), ]
		rownames(blk) <- NULL
		
		for (j in 1:nrow(crd)) {
			pos = crd$pos[j]
			
			crd$cov[j] = length(intersect(which(pos >= blk$l), which(pos <= blk$r)))
		}
		
		overl <- rbind(overl, crd)
	}
	
	overl
}


percent.cover <- function(x, pmin, pmax) {
	dist = pmax - pmin
	
	cut = cbind(x[order(x$pos), ], seg=NA)
	seg = 0
	
	for (i in 1:nrow(cut)) {
		if (cut$cov[i] == 0) {
			seg = seg + 1
			next
		}
		cut$seg[i] = seg
	}
	
	cut = split(cut, cut$seg)
	cut = lapply(cut, function(x) max(x$pos) - min(x$pos))
	
	covd = sum(unlist(cut))
	(covd / dist) * 100
}


table.coverage <- function(file) {
	
	cat(sprintf("%s\n", file))
	
	target <- sub(".+__(.+)__(.+)__(.+)\\..+", "\\1 \\2 \\3", file);
	target <- paste(rev(strsplit(target, " ", T)[[1]]), collapse = " ")
	
	tab <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
	
	tab <- cbind(tab,
							 sid = as.numeric(sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\1", tab$sample)),
							 key = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\2", tab$sample),
							 pop = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\3", tab$sample),
							 grp = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\4", tab$sample),
							 stringsAsFactors = FALSE)
	
	pmin = min(split.col(tab$lpos))
	pmax = max(split.col(tab$rpos))
	
	dat = overlap(tab)
	
	pcc = c()
	ihg = split(dat, dat$ihg)
	for (tag in sort(names(ihg))) {
		pcc = c(pcc, percent.cover(ihg[[tag]], pmin, pmax))
	}
	
	names(pcc) <- sprintf("cumulative.ihg.%d", 0:(length(pcc) - 1))
	
	cov = data.frame(sample = target, as.list(pcc), stringsAsFactors = FALSE)
	
	splt <- split(tab, tab$f)
	
	for (f in names(splt)) {
		dat = overlap(splt[[f]])
		
		pcc = c()
		ihg = split(dat, dat$ihg)
		for (tag in sort(names(ihg))) {
			pcc = c(pcc, percent.cover(ihg[[tag]], pmin, pmax))
		}
		
		names(pcc) <- sprintf("k.%s.ihg.%d", f, 0:(length(pcc) - 1))
		
		cov = cbind(cov, as.list(pcc), stringsAsFactors = FALSE)
	}
	
	cov
}


plot.coverage <- function(file) {
	
	cat(sprintf("%s\n", file))
	
	target <- sub(".+__(.+)__(.+)__(.+)\\..+", "\\1 \\2 \\3", file);
	target <- paste(rev(strsplit(target, " ", T)[[1]]), collapse = " ")
	
	tab <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
	
	tab <- cbind(tab,
							 sid = as.numeric(sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\1", tab$sample)),
							 key = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\2", tab$sample),
							 pop = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\3", tab$sample),
							 grp = sub("^ID-([0-9]+)__(.+)__(.+)__(.+)$", "\\4", tab$sample),
							 stringsAsFactors = FALSE)
	
	pmin = min(split.col(tab$lpos))
	pmax = max(split.col(tab$rpos))
	
	xseq <- seq(round(pmin/1000000)*1000000, round(pmax/1000000)*1000000, by = 1000000)
	xseq <- xseq[-1]
	
	
	string.cover <- function(x) {
		pcc = c()
		ihg = split(x, x$ihg)
		for (tag in sort(names(ihg))) {
			pcc = c(pcc, percent.cover(ihg[[tag]], pmin, pmax))
		}
		pcc = sprintf("IHG+%d = %0.3f %%", 0:(length(pcc) - 1), pcc)
		pcc = paste(pcc, collapse = ",  ")
	}
	
		
	dat = overlap(tab)
	cover <- cbind(k = sprintf(" Cumulative for k = {2, ..., %d}      Coverage:  %s", max(tab$f), string.cover(dat)), dat)
	
	splt <- split(tab, tab$f)
		
	for (f in names(splt)) {
		dat = overlap(splt[[f]])
		cover <- rbind(cover, cbind(k = sprintf(" k = %s      Coverage:  %s", f, string.cover(dat)), dat))
	}
	
	cover <- cbind(cover, shared = NA)
	i <- which(cover$pos %in% unique(tab$mpos))
	cover$shared[i] = cover$pos[i]
	
	
	yseq <- 0:max(cover$cov)
	yseq <- yseq[which(yseq / 10 == round(yseq / 10))]
	
	
	shared <- rbind(data.frame(x = unique(tab$mpos), y = -2 + 1), 
									data.frame(x = unique(tab$mpos), y = -2 - 1))
	
	
	ggplot(data = cover) + 
		theme_classic() + 
		theme(legend.position="none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(hjust = 0)) +
		
		facet_wrap(~k, ncol = 1) + 
		
		geom_hline(yintercept = yseq, colour = "grey95") +
		#geom_hline(yintercept = 1:length(unique(tab$f)), colour = "grey95", linetype="dashed") +
		
		geom_polygon(aes(x=pos, y=cov, group = ihg), fill="grey50", alpha = 1/3) + 
		
		geom_segment(aes(x = shared-0.5, xend=shared+0.5, y = -1, yend = -3), na.rm=TRUE, color="black", size=0.25, alpha = 1/3) +
		
		scale_x_continuous(expand = c(0, 0), breaks = xseq, labels = xseq / 1000000) +
		scale_y_continuous(expand = c(0, 0.5)) +
		xlab("Physical position (Mb)") + ylab("Regional sharing size, sub-sample sharing variant with target individual") + ggtitle(target)
	
	outfile <- sprintf("plot_coverage.%s.pdf", file)
	
	ggsave(filename = outfile, width = 25, height = 15, dpi = 300, limitsize = FALSE)
	
}



files <- dir(pattern = "^1kg.+\\.pairws$")[1:10]

coverage = NULL

for (file in files) {
	plot.pairws(file)
	plot.coverage(file)
	coverage = rbind(coverage, table.coverage(file))
}

save(coverage, file = "table.coverage.RData")

