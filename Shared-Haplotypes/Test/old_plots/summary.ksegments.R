


data <- read.table("./biallelic//summary.ksegments.txt", header = TRUE, stringsAsFactors = FALSE)


parse.kseg <- function(kseg)
{
	p <- strsplit(kseg, "|", TRUE)[[1]]
	p <- strsplit(p, ",|:")
	p <- lapply(p, as.numeric)
	
	pos <- sort(unique(c(unlist(lapply(p, function(x) x[2])), unlist(lapply(p, function(x) x[3])))))
	
	out <- data.frame(pos = pos, k=0)
	
	for (q in p)
	{
		min <- min(q[2:3])
		max <- max(q[2:3])
		
		k <- which(out$pos >= min & out$pos <= max)
		out$k[k] <- out$k[k] + q[1]
	}
	
	return(out)
}

parse.cseg <- function(cseg, root)
{
	p <- as.numeric(strsplit(cseg, ",|:")[[1]])
	
	min <- min(p[2:3])
	max <- max(p[2:3])
	
	return(data.frame(pos = c(min, root, max), k=p[1]))
}


splt <- split(data, data$n_samples)


library("parallel")

res <- mclapply(names(splt)[-(1:2)], function(f, splt) {
	
	s <- splt[[f]]
	
	cat(f, ":", nrow(s), "\n")
	
	out <- list()
	
	for (i in 1:nrow(s))
	{
		l <- r <- c()
		
		if (s$l_kseg[i] != ".")
			l <- parse.kseg(s$l_kseg[i])
		
		if (s$r_kseg[i] != ".")
			r <- parse.kseg(s$r_kseg[i])
		
		c <- parse.cseg(s$c_kseg[i], s$root_pos[i])
		
		out[[as.character(s$root_pos[i])]] <- rbind(l, c, r)
	}
	
	return(out)
	
} , splt, mc.cores=7)

names(res) <- names(splt)[-(1:2)]


save(res, file = "summary.ksegments.RData")



########


for (f in names(res))
{
	cat(f, "\n")
	
	for (p in names(res[[f]]))
	{
		pos <- as.numeric(p)
		
		c <- which(res[[f]][[p]]$pos == pos)
		centre <- res[[f]][[p]]$pos[c]
		
		if (length(c) > 1) {
			res[[f]][[p]] <- res[[f]][[p]][-(2:length(c)), ]
		}
		
		if (length(c) == 0) {
			res[[f]][[p]] <- NULL
			next
		}
		
		res[[f]][[p]]$pos <- res[[f]][[p]]$pos - centre[1]
	}
}



library("reshape2")
library("ggplot2")
library("scales")

for (nsam in c(100, 1000)) {
for (f in names(res))
{
	cat(f, "\n")
	
	r <- res[[f]]
	
	gg <- ggplot() + theme_bw() + theme_classic() + scale_y_continuous(expand = c(0,0)) + 
		xlab("Relative position") + ylab("k") + 
		scale_x_continuous(labels = comma) + 
		expand_limits(y=0) + 
		annotate("text",  x=Inf, y = Inf, label = sprintf("f%s", f), vjust=1, hjust=1)
	
	for (p in sample(names(r), nsam))
	{
		l_side <- which(r[[p]]$pos <= 0)
		r_side <- which(r[[p]]$pos >= 0)
		
		if (length(l_side) > 0)
			gg <- gg + geom_step(data = r[[p]][l_side, ], aes(x=pos, y=k), alpha = 0.05, direction = "hv") + geom_point(data=r[[p]][1, ], mapping=aes(x=pos, y=k), size=2, alpha = 0.05)
		
		if (length(r_side) > 0)
			gg <- gg + geom_step(data = r[[p]][r_side, ], aes(x=pos, y=k), alpha = 0.05, direction = "vh") + geom_point(data=r[[p]][nrow(r[[p]]), ], mapping=aes(x=pos, y=k), size=2, alpha = 0.05)
	}
	
	ggsave(filename = sprintf("summary.ksegments.biallelic.plot%d.f%s.png", nsam, f), width=15, height=5)
}
}









