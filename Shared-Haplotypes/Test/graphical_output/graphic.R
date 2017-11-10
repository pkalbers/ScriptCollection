#
# graphical output of shared haplotype structures
#

library("ggplot2")


raw <- read.table("struct.print.txt", header = TRUE, stringsAsFactors = FALSE)
mrk <- read.table("struct.marker", header = TRUE, stringsAsFactors = FALSE)

foc <- sample(which(raw$site_count == 5 & raw$sample_n == 5), 1)
shs <- raw[foc, ]


bound.pos.min <- min(mrk$position)
bound.pos.max <- max(mrk$position)
bound.map.min <- min(mrk$map_dist)
bound.map.max <- max(mrk$map_dist)


parse.comma <- function(r)
{
	unlist(strsplit(r, split = ",", fixed = TRUE))
}


shs.core.focus <- function(shs)
{
	list(hap = shs$site_hap,
			 rac = shs$site_count,
			 raf = shs$site_freq,
			 n_sample = shs$sample_n,
			 c_pos = shs$site_pos,
			 c_map = shs$site_map,
			 l_pos = shs$core_pos_from,
			 l_map = shs$core_map_from,
			 r_pos = shs$core_pos_to,
			 r_map = shs$core_map_to)
}

shs.core.sample <- function(shs)
{
	l <- unlist(strsplit(shs$break_gen_l, ":", TRUE))
	r <- unlist(strsplit(shs$break_gen_r, ":", TRUE))
	
	data.frame(id  = parse.comma(shs$sample_id),
						 key = parse.comma(shs$sample_key),
						 pop = parse.comma(shs$sample_pop),
						 grp = parse.comma(shs$sample_grp),
						 c_gen = parse.comma(shs$site_gen),
						 l_gen = parse.comma(l[3]),
						 r_gen = parse.comma(r[3]),
						 stringsAsFactors = FALSE)
}

shs.core.break <- function(shs)
{
	l <- unlist(strsplit(shs$break_gen_l, ":", TRUE))
	r <- unlist(strsplit(shs$break_gen_r, ":", TRUE))
	
	list(l_pos = as.numeric(l[1]),
			 l_map = as.numeric(l[2]),
			 r_pos = as.numeric(r[1]),
			 r_map = as.numeric(r[2]))
}

shs.core <- function(shs)
{
	list(focus = shs.core.focus(shs), sample = shs.core.sample(shs), breaks = shs.core.break(shs))
}


shs.tree <- function(shs)
{
	shs.lr.tree <- function(shs, tree)
	{
		t <- unlist(strsplit(tree, "|", TRUE))
		t <- strsplit(t, ":", TRUE)
		
		d <- data.frame()
		for (i in 1:length(t))
		{
			loc <- unlist(strsplit(t[[i]][2], "_", TRUE))
			
			pos.from  <- as.numeric(t[[i]][3])
			pos.to    <- as.numeric(t[[i]][4])
			pos.break <- if (t[[i]][9] == ".") { if (pos.from > pos.to) bound.pos.min else bound.pos.max } else { as.numeric(t[[i]][9]) }
			
			map.from  <- as.numeric(t[[i]][5])
			map.to    <- as.numeric(t[[i]][6])
			map.break <- if (t[[i]][9] == ".") { if (map.from > map.to) bound.map.min else bound.map.max } else { as.numeric(t[[i]][10]) }
			
			d <- rbind(d, data.frame(connect_to = t[[i]][1],
															 connect_as = loc[1],
															 indic = as.numeric(loc[2]),
															 level = as.numeric(loc[3]),
															 pos_from = pos.from,
															 pos_to = pos.to,
															 pos_break = pos.break,
															 map_from = map.from,
															 map_to = map.to,
															 map_break = map.break,
															 n_subsample = as.numeric(t[[i]][7]), 
															 stringsAsFactors = FALSE))
		}
		
		sample <- shs.core.sample(shs)
		sample$c_gen <- NULL
		sample$l_gen <- NULL
		sample$r_gen <- NULL
		
		e <- list()
		for (i in 1:length(t))
		{
			s <- data.frame()
			
			ids <- as.numeric(parse.comma(t[[i]][8]))
			
			for (k in 1:length(ids))
			{
				x <- which(sample$id == ids[k])
				
				if (length(x) == 0)
					stop("Sample ID not found")
				
				s <- rbind(s, sample[x, ])
			}
			
			if (t[[i]][9] == ".")
				s <- cbind(s, gen = "-/-")
			else
				s <- cbind(s, gen = parse.comma(t[[i]][11]))
			
			e[[as.character(d$connect_as[i])]] <- s
		}
		
		list(tree = d, subsample = e)
	}
	
	l <- r <- NA
	
	if (shs$l_tree != ".")
		l <- shs.lr.tree(shs, shs$l_tree)
	
	if (shs$l_tree != ".")
		r <- shs.lr.tree(shs, shs$r_tree)
	
	list(l = l, r = r)
}



sample.colors <- function(shs)
{
	key <- parse.comma(shs$sample_key)
	pop <- parse.comma(shs$sample_pop)
	grp <- parse.comma(shs$sample_grp)
	let <- sample.letter(shs)
	col <- c()
	
	for (i in 1:length(key))
		col <- c(col, paste(sprintf("%s:", let[(key[i])]), key[i], pop[i], grp[i]))

	names(col) <- key
	
	col
}

sample.letter <- function(shs)
{
	key <- parse.comma(shs$sample_key)
	let <- LETTERS[1:(length(key))]
	names(let) <- rev(key)
	
	let
}



yloc <- function(n, center = 0)
{
	y <- 1:n - mean(1:n) + center
	p <- c(min(y) - 1, max(y) + 1)
	
	list(y.coord = y, padding = p, center = center)
}


tree.coord <- function(shs)
{
	core <- shs.core(shs)
	tree <- shs.tree(shs)
	
	enum.tips <- function(t)
	{
		tips <- data.frame()
		
		ends <- which(! t$connect_as %in% t$connect_to)
		
		for (i in ends)
		{
			r <- t$connect_to[i]
			n <- length(which(t$connect_to == r))
			
			x <- 2
			
			#if (n == 1) {
				if (r == "ROOT") {
					x <- core$focus$n_sample
				} else {
					x <- t$n_subsample[which(t$connect_as == r)]
				}
			#}
			
			tips <- rbind(tips, data.frame(tag = t$connect_as[i], n = t$n_subsample[i], x = x, stringsAsFactors = FALSE))
		}
		tips
	}
	
	l <- enum.tips(tree$l$tree)
	r <- enum.tips(tree$r$tree)
	
	l <- cbind(l, npad = l$x + 2, center = NA)
	r <- cbind(r, npad = r$x + 2, center = NA)
	
	last <- mean(1:(sum(l$npad)))
	for (i in 1:nrow(l)) {
		l$center[i] <- last - mean(1:(l$npad[i]))
		last <- last - l$npad[i]
	}
	
	last <- mean(1:(sum(r$npad)))
	for (i in 1:nrow(r)) {
		r$center[i] <- last - mean(1:(r$npad[i]))
		last <- last - r$npad[i]
	}
	
	l.coord <- list()
	for (tag in l$tag) {
		i <- which(l$tag == tag)
		l.coord[[tag]] <- yloc(l$n[i], l$center[i])
	}
	
	r.coord <- list()
	for (tag in r$tag) {
		i <- which(r$tag == tag)
		r.coord[[tag]] <- yloc(r$n[i], r$center[i])
	}
	
	sub.coord <- function(t, coord)
	{
		while(length(coord) != nrow(t))
		{
			for (tag in names(coord))
			{
				r <- t$connect_to[which(t$connect_as == tag)]
				
				if (r == "ROOT") next
				if (r %in% names(coord)) next
				
				i <- which(t$connect_to == r)
				
				if (all(t$connect_as[i] %in% names(coord)))
				{
					subs <- coord[t$connect_as[i]]
					center <- if (length(i) == 1) subs[[1]]$center else mean(range(unlist(lapply(subs, function(x) x$y.coord))))
					coord[[r]] <- yloc(t$n_subsample[which(t$connect_as == r)], center)
				}
				
			}
		}
		
		coord
	}
	
	l.coord <- sub.coord(tree$l$tree, l.coord)
	r.coord <- sub.coord(tree$r$tree, r.coord)
	
	list(l = l.coord, r = r.coord)
}




plot.core <- function(shs)
{
	core <- shs.core(shs)
	tree <- shs.tree(shs)
	keys <- parse.comma(shs$sample_key)
	cols <- sample.colors(shs)
	lets <- sample.letter(shs)
	
	ylocs <- yloc(core$focus$n_sample)
	
	thick.l <- sum(tree$l$tree[which(! tree$l$tree$connect_as %in% tree$l$tree$connect_to), ]$n_subsample + 2)
	thick.r <- sum(tree$r$tree[which(! tree$r$tree$connect_as %in% tree$r$tree$connect_to), ]$n_subsample + 2)
	thick   <- min(c(3, 300 / max(c(thick.l, thick.r))))
	
	# break lines
	ld <- data.frame()
	for (i in 1:core$focus$n_sample)
	{
		ld <- rbind(ld, data.frame(x = c(core$breaks$l_pos / 1000, core$breaks$r_pos / 1000), y = ylocs$y.coord[i]))
	}
	ld <- geom_line(data = ld, mapping = aes(x=x, y=y, group=y), size = thick / 3, colour = "grey")
	
	# core lines
	ls <- data.frame()
	for (i in 1:core$focus$n_sample)
	{
		ls <- rbind(ls, data.frame(x = c(core$focus$l_pos / 1000, core$focus$r_pos / 1000), y = ylocs$y.coord[i], col = as.character(cols[keys[i]]), lab = as.character(lets[keys[i]])))
	}
	ls <- geom_line(data = ls, mapping = aes(x=x, y=y, group=y, colour = col), size = thick)
	
	# vertical line at position
	vl <- geom_vline(xintercept = core$focus$c_pos / 1000, size = thick / 3, linetype = "dotted")
	
	# focus points
	p <- data.frame(x = core$focus$c_pos / 1000, y = ylocs$y.coord)
	p <- geom_point(data = p, mapping = aes(x=x, y=y), shape=15, color="white", size = thick * 3, alpha=0.75)
	
	# annotation
	a <- annotate(geom = "text", x = core$focus$c_pos / 1000, y = ylocs$y.coord, label = lets[core$sample$key], hjust = 0.5, size = thick * 2)
	
	return(list(ld, ls, vl, p, a))
}


plot.tree <- function(shs)
{
	core <- shs.core(shs)
	tree <- shs.tree(shs)
	coord <- tree.coord(shs)
	cols <- sample.colors(shs)
	lets <- sample.letter(shs)

	thick.l <- sum(tree$l$tree[which(! tree$l$tree$connect_as %in% tree$l$tree$connect_to), ]$n_subsample + 2)
	thick.r <- sum(tree$r$tree[which(! tree$r$tree$connect_as %in% tree$r$tree$connect_to), ]$n_subsample + 2)
	thick   <- min(c(3, 300 / max(c(thick.l, thick.r))))
	
	plot.ublock <- function(t, s, c)
	{
		block <- data.frame()
		for (i in 1:nrow(t))
		{
			tag <- t$connect_as[i]
			loc <- c[[tag]]
			
			prv <- t$connect_to[i]
			prev.brk <- 0
			if (prv == "ROOT") {
				l <- core$breaks$l_pos
				r <- core$breaks$r_pos
				prev.brk <- if (as.numeric(dist(c(l, t$pos_break[i]))) < as.numeric(dist(c(r, t$pos_break[i])))) l else r
			} else {
				j <- which(t$connect_as == prv)
				prev.brk <- t$pos_break[j]
			}
			
			for (k in 1:(t$n_subsample[i])) 
				block <- rbind(block, data.frame(x = c(t$pos_break[i] / 1000, prev.brk / 1000), y = loc$y.coord[k], group = sprintf("%d_%d", i, k)))
		}
		
		geom_line(data = block, mapping = aes(x=x, y=y, group=group), size = thick / 3, colour = "grey")
	}
	
	plot.block <- function(t, s, c)
	{
		block <- list()
		for (i in 1:nrow(t)) 
		{
			tag <- t$connect_as[i]
			loc <- c[[tag]]
			
			for (k in 1:(t$n_subsample[i]))
				block <- rbind(block, data.frame(x = c(t$pos_from[i] / 1000, t$pos_to[i] / 1000), y = loc$y.coord[k], group = sprintf("%d_%d", i, k), col = as.character(cols[s[[tag]]$key[k]])))
		}
		geom_line(data = block, mapping = aes(x=x, y=y, group=group, colour = col), size = thick)
	}
	
	plot.break <- function(t, s, c)
	{
		bound <- -1
		
		tips <- t[which(! t$connect_as %in% t$connect_to), ]$connect_as
		
		block <- data.frame()
		genbr <- data.frame()
		
		for (i in 1:nrow(t))
		{
			tag <- t$connect_as[i]
			loc <- c[[tag]]
			
			prev <- t$connect_to[i]
			
			if (prev %in% block$tag)
				next
			
			prev.brk <- 0
			if (prev == "ROOT") {
				l <- core$breaks$l_pos
				r <- core$breaks$r_pos
				prev.brk <- if (as.numeric(dist(c(l, t$pos_break[i]))) < as.numeric(dist(c(r, t$pos_break[i])))) l else r
			} else {
				prev.brk <- t$pos_break[which(t$connect_as == prev)]
			}
			
			k <- which(t$connect_to == t$connect_to[i])
			shrd <- t$connect_as[k]
			prev.y <- unlist(lapply(c[shrd], function(x) x$padding))
			
			if (prev == "ROOT")
				prev.y <- c(prev.y, yloc(core$focus$n_sample)$padding)
			
			prev.y <- range(prev.y)
			
			block <- rbind(block, data.frame(x = prev.brk / 1000, y = prev.y, tag = prev))
			
			
			if (prev == "ROOT") {
				prev.y <- yloc(core$focus$n_sample)$y.coord
				l <- core$breaks$l_pos
				r <- core$breaks$r_pos
				prev.g <- if (as.numeric(dist(c(l, t$pos_break[i]))) < as.numeric(dist(c(r, t$pos_break[i])))) core$sample$l_gen else core$sample$r_gen
				
			} else {
				prev.y <- c[[prev]]$y.coord
				prev.g <- s[[prev]]$gen
			}
			
			prev.g <- as.character(unlist(lapply(strsplit(as.character(prev.g), "/", TRUE), function(x) sum(as.numeric(x)) )))
			
			genbr <- rbind(genbr, data.frame(x = prev.brk / 1000, y = prev.y, gen = prev.g))

			if (t$pos_break[i] == bound.pos.min)
				bound <- bound.pos.min
			if (t$pos_break[i] == bound.pos.max)
				bound <- bound.pos.max
		}
		
		bound.l <- bound.r <- NULL
		
		if (bound == bound.pos.min)
			bound.l <- geom_vline(xintercept = bound.pos.min / 1000, size = thick / 3, colour = "grey", linetype = "dashed")
		if (bound == bound.pos.max)
			bound.r <- geom_vline(xintercept = bound.pos.max / 1000, size = thick / 3, colour = "grey", linetype = "dashed")
		
		list(geom_line(data = block, mapping = aes(x=x, y=y, group=tag), size = thick / 3, colour = "grey"), bound.l, bound.r,
				 geom_text(data = genbr, mapping = aes(x=x, y=y, label = gen), size = thick * 2, fontface="bold", colour = "white"),
				 geom_text(data = genbr, mapping = aes(x=x, y=y, label = gen), size = thick * 2)
				 )
	}
	
	plot.annot <- function(t, s, c)
	{
		annot <- data.frame()
		gentp <- data.frame()
		
		tip <- t[which(! t$connect_as %in% t$connect_to), ]
		
		right <- if (tip$pos_break[1] > tip$pos_to[1]) FALSE else TRUE
		
		for (i in 1:nrow(tip))
		{
			tag <- tip$connect_as[i]
			keys <- s[[tag]]$key
			
			annot <- rbind(annot, data.frame(x = tip$pos_break[i] / 1000, y = c[[tag]]$y.coord, label = lets[keys]))
			gentp <- rbind(gentp, data.frame(x = tip$pos_break[i] / 1000, y = c[[tag]]$y.coord, label = as.character(unlist(lapply(strsplit(as.character(s[[tag]]$gen), "/", TRUE), function(x) sum(as.numeric(x)) ))) ))
		}
		
		list(annotate(geom = "text", x = annot$x, y = annot$y, label = annot$label, hjust = if (right) 1.5 else -0.5, size = thick * 2), 
				 annotate(geom = "text", x = gentp$x, y = gentp$y, label = gentp$label, size = thick * 2, fontface="bold", colour = "white"),
				 annotate(geom = "text", x = gentp$x, y = gentp$y, label = gentp$label, size = thick * 2))
	}
	
	c(plot.ublock(tree$l$tree, tree$l$subsample, coord$l), 
		plot.ublock(tree$r$tree, tree$r$subsample, coord$r),
		plot.block(tree$l$tree, tree$l$subsample, coord$l), 
		plot.block(tree$r$tree, tree$r$subsample, coord$r),
		plot.break(tree$l$tree, tree$l$subsample, coord$l), 
		plot.break(tree$r$tree, tree$r$subsample, coord$r),
		plot.annot(tree$l$tree, tree$l$subsample, coord$l),
		plot.annot(tree$r$tree, tree$r$subsample, coord$r)
		)
}

plot.legend <- function(shs)
{
	tree <- shs.tree(shs)
	coord <- tree.coord(shs)
	cols <- sample.colors(shs)
	lets <- sample.letter(shs)

	spl <- data.frame(key = parse.comma(shs$sample_key),
										pop = parse.comma(shs$sample_pop),
										grp = parse.comma(shs$sample_grp),
										col = NA,
										let = NA,
										stringsAsFactors = FALSE
										)
	
	for (i in 1:nrow(spl)) {
		spl$let[i] <- lets[(spl$key[i])]
		spl$col[i] <- cols[(spl$key[i])]
	}
	
	spl <- spl[order(spl$let), ]
	
	thick.l <- sum(tree$l$tree[which(! tree$l$tree$connect_as %in% tree$l$tree$connect_to), ]$n_subsample + 2)
	thick.r <- sum(tree$r$tree[which(! tree$r$tree$connect_as %in% tree$r$tree$connect_to), ]$n_subsample + 2)
	thick   <- min(c(3, 300 / max(c(thick.l, thick.r))))
	right   <- (thick.l > thick.r)
	
	y <- if (right) max(unlist(lapply(coord$l, function(x) x$padding))) else max(unlist(lapply(coord$r, function(x) x$padding)))
	
	spl <- cbind(spl, 
							 x = shs$site_pos / 1000, 
							 y = y - ((0:(nrow(spl) - 1))),
							 label = paste(spl$let, "-", spl$key, spl$pop, spl$grp)
							 )
	
	geom_text(data = spl, mapping = aes(x=x, y=y, label=label), alpha = 0.5, family="mono")
}


plot.f <- function(f)
{
	foc <- sample(which(raw$site_count == f & raw$sample_n == f), 1)
	shs <- raw[foc, ]
	
	g <- ggplot() + theme_classic() + theme(axis.title.y = element_blank(), 
																					axis.text.y = element_blank(), 
																					axis.ticks.y = element_blank(), 
																					axis.line.y = element_blank(),
																					axis.title.x = element_text(size = 25), 
																					axis.text.x = element_text(size = 25), 
																					axis.ticks.x = element_line(size = 2), 
																					axis.line.x = element_line(size = 2), 
																					plot.title = element_text(size = 35, face="bold"),
																					legend.title=element_blank()
	) + xlab("Position [Kb]") + ggtitle(sprintf("f %d", shs$sample_n))
	
	g <- g + plot.core(shs) + plot.tree(shs)
	
	ggsave(plot = g, filename = sprintf("plot.f%d.%d.pdf", shs$sample_n, shs$site_pos), width = 50, height = 35, limitsize=FALSE)
}

for (f in c(5, 10, 15, 20)) {
	plot.f(f)
	plot.f(f)
	plot.f(f)
	plot.f(f)
	plot.f(f)
}





