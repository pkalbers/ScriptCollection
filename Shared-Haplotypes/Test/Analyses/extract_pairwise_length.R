

mrk <- read.table("struct.marker", header = TRUE, stringsAsFactors = FALSE)

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


shs.tree <- function(shs,
										 bound.pos.min, bound.pos.max,
										 bound.map.min, bound.map.max)
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






pair.dist <- function(shs,
											bound.pos.min, bound.pos.max,
											bound.map.min, bound.map.max) {
	
	if (shs$break_gen_l == "." || shs$break_gen_r == ".")
		return(NULL)
	
	core <- shs.core(shs)
	tree <- shs.tree(shs, 
									 bound.pos.min, bound.pos.max,
									 bound.map.min, bound.map.max)
	
	if(length(tree$l$subsample) == 0)
		return(NULL)
		
	indv <- Reduce(rbind, tree$l$subsample)
	indv$gen <-NULL
	indv <- unique(indv)
	rownames(indv) <- indv$id
	
	
	node.l <- list()
	node.r <- list()
	for (s in indv$id) {
		node <- c()
		for (tag in names(tree$l$subsample)) {
			if (s %in% tree$l$subsample[[tag]]$id) node <- c(node, tag)
		}
		node.l[[s]] <- node
		
		node <- c()
		for (tag in names(tree$r$subsample)) {
			if (s %in% tree$r$subsample[[tag]]$id) node <- c(node, tag)
		}
		node.r[[s]] <- node
	}
	
	comb <- combn(indv$id, 2)
	
	pair <- data.frame(a = comb[1, ], b = comb[2, ], 
										 lpos = NA, rpos = NA, 
										 lmap = NA, rmap = NA, 
										 a.pop = NA, b.pop = NA,
										 a.grp = NA, b.grp = NA,
										 stringsAsFactors = FALSE)
	
	for (i in 1:nrow(pair)) {
		
		lpos <- rpos <- c()
		lmap <- rmap <- c()
		
		nl <- intersect(node.l[[ pair$a[i] ]], node.l[[ pair$b[i] ]])
		nr <- intersect(node.r[[ pair$a[i] ]], node.r[[ pair$b[i] ]])
		
		if (length(nl) == 0) {
			lpos <- c(lpos, core$focus$l_pos)
			lmap <- c(lmap, core$focus$l_map)
		} else {
			for (con in nl) {
				tl <- tree$l$tree[which(tree$l$tree$connect_as == con), ]
				lpos <- c(lpos, tl$pos_to)
				lmap <- c(lmap, tl$map_to)
			}
		}
		
		if (length(nr) == 0) {
			rpos <- c(rpos, core$focus$r_pos)
			rmap <- c(rmap, core$focus$r_map)
		} else {
			for (con in nr) {
				tr <- tree$r$tree[which(tree$r$tree$connect_as == con), ]
				rpos <- c(rpos, tr$pos_to)
				rmap <- c(rmap, tr$map_to)
			}
		}
		
		if (length(lpos) == 0 || length(rpos) == 0)
			return(NULL)
		
		pair$lpos[i] <- min(lpos)
		pair$rpos[i] <- max(rpos)
		pair$lmap[i] <- min(lmap)
		pair$rmap[i] <- max(rmap)
	}
	
	pair$a.pop <- indv[pair$a, "pop"]
	pair$b.pop <- indv[pair$b, "pop"]
	
	pair$a.grp <- indv[pair$a, "grp"]
	pair$b.grp <- indv[pair$b, "grp"]
	
	return(pair)
}

library("parallel")



raw <- read.table("struct.print.txt", header = TRUE, stringsAsFactors = FALSE)



k = 5

focal <- which(raw$site_count == k & raw$sample_n == k)
print(length(focal))

shs <- list()
for (foc in focal)
	shs[[ as.character(foc) ]] <- raw[foc, ]

save(shs, file = "_tmp.shs05.RData")



k = 25

focal <- which(raw$site_count == k & raw$sample_n == k)
print(length(focal))

shs <- list()
for (foc in focal)
	shs[[ as.character(foc) ]] <- raw[foc, ]

save(shs, file = "_tmp.shs25.RData")


rm(raw, shs)


load("_tmp.shs05.RData")
res <- mclapply(shs, pair.dist, bound.pos.min, bound.pos.max, bound.map.min, bound.map.max, mc.cores = 30)
res <- Reduce(rbind, res)
save(res, file = "_pairwise.05.RData")

rm(shs)

load("_tmp.shs25.RData")
res <- mclapply(shs, pair.dist, bound.pos.min, bound.pos.max, bound.map.min, bound.map.max, mc.cores = 30)
res <- Reduce(rbind, res)
save(res, file = "_pairwise.05.RData")

#median(res$rpos - res$lpos)
#median(res$rmap - res$lmap)


stop()




load("_pairwise.05.RData")
load("_pairwise.25.RData")


length(which(res$b.grp != res$a.grp)) / nrow(res)


sam <- unique(c(res$a, res$b))
names(sam) <- sam
splt <- lapply(sam, function(x, r) r[ unique(c(which(r$a == x), which(r$b == x))) , ] , res )
z <- lapply(splt, function(x) length(which(x$a.grp != x$b.grp)) / nrow(x) )
median(unlist(z))


pop <- data.frame(sam = sam, pop = NA, stringsAsFactors = FALSE)
for (i in 1:nrow(pop)) {
	pop$pop[i] <- splt[[ pop$sam[i] ]]$a.grp[1]
}
sub <- split(pop, pop$pop)
lapply(sub, function(x, s) { y <- s[x$sam]; median(unlist(lapply(y, function(x) length(which(x$a.grp != x$b.grp)) / nrow(x) ))) }, splt)



splt <- split(res, res$b)
z <- lapply(splt, function(x) length(which(x$a.grp != x$b.grp)) / nrow(x) )
median(unlist(z))



res <- cbind(res, pair = paste(res$a, res$b))

splt <- split(res, res$pair)




