
args <- commandArgs(TRUE)

read.table(args[1], TRUE)->stat
read.table(args[2], TRUE)->shap


library("parallel")


accum.incons <- function(i, coor, stat, max, intv=150) {
	a <- coor[i, "from"]
	o <- coor[i, "mid"]
	z <- coor[i, "to"]
	
	# excluding haplotypes at chromosome boundaries
	if (stat$position[1] == a || stat$position[nrow(stat)] == z)
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
	
	# pyhsical positions
	#phypos.l <- stat$position[l.in]
	#phypos.p <- stat$position[io]
	#phypos.r <- stat$position[r.in]
	
	# get inconsistency prob.
	incs.l <- ( stat$n_inconsistencies[l] * stat$maf[l] ) / max / length(l)
	incs.r <- ( stat$n_inconsistencies[r] * stat$maf[r] ) / max / length(r)
	
	# genetic distances
	gendist.l <- l - io
	gendist.r <- r - io
	gendist.l <- gendist.l / -min(gendist.l)
	gendist.r <- gendist.r /  max(gendist.r)
	
	# physical distances
	#phydist.l <- phypos.l - phypos.p
	#phydist.r <- phypos.r - phypos.p
	#phydist.l <- phydist.l / -min(phydist.l)
	#phydist.r <- phydist.r /  max(phydist.r)
	
	cut.l <- cut(gendist.l, seq(0, -1, length.out=intv+1), include.lowest=TRUE, right=FALSE)
	incs.l <- split(incs.l, cut.l)
	incs.l <- lapply(incs.l, sum)#, na.rm=TRUE)
	
	cut.r <- cut(gendist.r, seq(0, 1, length.out=intv+1), include.lowest=TRUE, right=TRUE)
	incs.r <- split(incs.r, cut.r)
	incs.r <- lapply(incs.r, sum)#, na.rm=TRUE)
	
	return(c(unlist(incs.l), unlist(incs.r)))
}



coor <- cbind(from=shap$from, mid=shap$position, to=shap$to)

n <- nrow(coor)

max <- choose(2657, 2) * 0.5

full <- 1:n
chunks <- split(full, ceiling(seq_along(full)/100000))

map <- c()

i <- 0
for (chunk in chunks) {	
	i <- i + 1
	cat(i, "/", length(chunks), "\n")
	
	result <- mclapply(chunk, accum.incons, coor, stat, max, mc.cores=3)
	
	col <- max(unlist(lapply(result, length)))
	mat <- matrix(unlist(result), ncol=col, byrow=TRUE)
	map <- rbind(map, mat)
	
	save(map, file=sprintf("map.%s.RData", basename(args[1])))
}


stop("DONE!")


est <- apply(map, 2, mean, na.rm=TRUE)


w <- length(est)
x <- seq(-1.5, 1.5, length.out=w)

est[1] <- NA
est[w] <- NA
est[which(x > -0.01 & x < 0.01)] <- NA

plot(x, est, type='l', lwd=1, 
	 xlab="Relative distance to rare variant", ylab="Mean aggregated inconsistency")
x.c <- x[which(x > -0.01 & x < 0.01)]
x.l <- x[which(x < -0.99 & x > -1.01)]
x.r <- x[which(x > 0.99 & x < 1.01)]
polygon(c(x.l, rev(x.l)), c(1,1,-1,-1), border=NA, col="grey")
polygon(c(x.r, rev(x.r)), c(1,1,-1,-1), border=NA, col="grey")
polygon(c(x.c, rev(x.c)), c(1,1,-1,-1), border=NA, col="grey")
lines(x, est, lwd=2)
box(lwd=2)




#cin <- apply(map, 2, function(x) (sd(x, na.rm=TRUE) / sqrt(length(na.omit(x)))) * 1.96)
#cin[(w/2):((w/2)+1)] <- NA
#cin[1] <- NA
#cin[w] <- NA
# x.l <- x[2:((w/2)-1)]
# x.r <- x[((w/2)+2):(w-1)]
# est.l <- est[2:((w/2)-1)]
# est.r <- est[((w/2)+2):(w-1)]
# cin.l <- cin[2:((w/2)-1)]
# cin.r <- cin[((w/2)+2):(w-1)]
# polygon(c(x.l, rev(x.l)), c(est.l+cin.l, rev(est.l-cin.l)), border=NA, col="grey")
# polygon(c(x.r, rev(x.r)), c(est.r+cin.r, rev(est.r-cin.r)), border=NA, col="grey")
#lines(x, est, lwd=1.25)








