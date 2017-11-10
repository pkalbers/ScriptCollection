
args <- commandArgs(TRUE)

read.table(args[1], TRUE)->stat
read.table(args[2], TRUE)->shap


library("parallel")


accum.incons <- function(i, coor, stat, intv=150) {
	a <- coor$from[i]
	o <- coor$mid[i]
	z <- coor$to[i]
	
	# excluding haplotypes at chromosome boundaries
	if (min(stat$position) == a || max(stat$position) == z)
		return(NULL)
	
	# getting vector indices
	ia <- which(stat$position == a)
	io <- which(stat$position == o)
	iz <- which(stat$position == z)
	
	# left and right from shared variant
	l <- (io-1):ia
	r <- (io+1):iz
	l <- c(l, (ia-1):(ia - floor(length(l) / 2)))
	r <- c(r, (iz+1):(iz + floor(length(r) / 2)))
	
	if (min(l) < 1 || max(r) > nrow(stat))
		return(NULL)
	
	if (length(l) < intv || length(r) < intv)
		return(NULL)
	
	# get data
	inc.l <- stat$n_inconsistencies[l] * stat$maf[l]
	inc.r <- stat$n_inconsistencies[r] * stat$maf[r]
	
	# physical positions
	phypos <- stat$position[io]
	phypos.l <- stat$position[l]
	phypos.r <- stat$position[r]
	
	# physical distances
	phydist.l <- phypos.l - phypos
	phydist.r <- phypos.r - phypos
	phydist.l <- phydist.l / -min(phydist.l)
	phydist.r <- phydist.r /  max(phydist.r)
	
	cut.l <- cut(phydist.l, seq(0, -1, length.out=intv+1), include.lowest=TRUE, right=FALSE)
	res.l <- lapply(split(inc.l, cut.l), sum)
	
	cut.r <- cut(phydist.r, seq(0, 1, length.out=intv+1), include.lowest=TRUE, right=TRUE)
	res.r <- lapply(split(inc.r, cut.r), sum)
	
	return(c(unlist(res.l), unlist(res.r)))
}



coor <- data.frame(from=shap$from, mid=shap$position, to=shap$to)

n <- nrow(coor)

#max <- choose(2657, 2) * 0.5

full <- 1:n
chunks <- split(full, ceiling(seq_along(full)/100000))

map <- c()

i <- 0
for (chunk in chunks) {	
	i <- i + 1
	cat(i, "/", length(chunks), "\n")
	
	result <- mclapply(chunk, accum.incons, coor, stat, mc.cores=24)
	
	col <- max(unlist(lapply(result, length)))
	mat <- matrix(unlist(result), ncol=col, byrow=TRUE)
	map <- rbind(map, mat)
	
	#save(map, file=sprintf("map.incons.%s.RData", basename(args[1])))
}
save(map, file=sprintf("map.incons.%s.RData", basename(args[1])))

stop("DONE!")


est <- apply(map, 2, mean, na.rm=TRUE)


w <- length(est)
x <- seq(-1.5, 1.5, length.out=w)

est[1] <- NA
est[w] <- NA
est[which(x > -0.01 & x < 0.01)] <- NA

plot(x, est, type='l', lwd=1, 
	 xlab="Relative distance to rare variant", ylab="Mean aggregated inconsistency rate")
x.c <- x[which(x > -0.01 & x < 0.01)]
x.l <- x[which(x < -0.99 & x > -1.01)]
x.r <- x[which(x > 0.99 & x < 1.01)]
polygon(c(x.l, rev(x.l)), c(1,1,-1,-1), border=NA, col="grey")
polygon(c(x.r, rev(x.r)), c(1,1,-1,-1), border=NA, col="grey")
polygon(c(x.c, rev(x.c)), c(1,1,-1,-1), border=NA, col="grey")
lines(x, est, lwd=2)
box(lwd=2)







