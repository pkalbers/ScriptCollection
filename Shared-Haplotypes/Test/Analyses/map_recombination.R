
args <- commandArgs(TRUE)

read.table(args[1], TRUE)->stat
read.table(args[2], TRUE)->shap
read.table(args[3], TRUE)->rec

rec$Rate.cM.Mb. <- runif(length(rec$Rate.cM.Mb.), min=0, max=1)

library("parallel")


accum.rec <- function (i, shap, stat, rec, intv=150) {
	a <- shap$from[i]
	o <- shap$position[i]
	z <- shap$to[i]
	
	# getting STAT indices
	ia <- which(stat$position == a)
	io <- which(stat$position == o)
	iz <- which(stat$position == z)
	
	# left and right from shared variant
	l <- (io):ia
	r <- (io):iz
	l <- c(l, (ia-1):(ia - floor(length(l) / 2)))
	r <- c(r, (iz+1):(iz + floor(length(r) / 2)))
	
	if (min(l) < 1 || max(r) > nrow(stat))
		return(NULL)
	
	if (length(l) < intv || length(r) < intv)
		return(NULL)
	
	l <- rev(l) # !!!
	
	# positions in STAT
	posl <- stat$position[l]
	posr <- stat$position[r]
	
	# getting REC indices
	il <- which(rec$Position.bp. %in% posl)
	ir <- which(rec$Position.bp. %in% posr)
	jl <- which(posl %in% rec$Position.bp.)
	jr <- which(posr %in% rec$Position.bp.)
	
	if (length(il) == 0 || length(ir) == 0)
		return(NULL)
	
	# get rates for region
	rate.l <- rep(0, length(posl))
	rate.r <- rep(0, length(posr))
	rate.l[jl] <- rec$Rate.cM.Mb.[il]
	rate.r[jr] <- rec$Rate.cM.Mb.[ir]
	
	# genetic distances
	gendist.l <- l - io
	gendist.r <- r - io
	gendist.l <- gendist.l / -min(gendist.l)
	gendist.r <- gendist.r /  max(gendist.r)
	
	cut.l <- cut(gendist.l, seq(0, -1, length.out=intv+1), include.lowest=TRUE, right=FALSE)
	rate.l <- split(rate.l, cut.l)
	rate.l <- lapply(rate.l, mean, na.rm=TRUE)
	
	cut.r <- cut(gendist.r, seq(0, 1, length.out=intv+1), include.lowest=TRUE, right=TRUE)
	rate.r <- split(rate.r, cut.r)
	rate.r <- lapply(rate.r, mean, na.rm=TRUE)
	
	return(c(unlist(rate.l), unlist(rate.r)))
}



n <- nrow(shap)
full <- 1:n
chunks <- split(full, ceiling(seq_along(full)/100000))

recmap <- c()

i <- 0
for (chunk in chunks) {	
	i <- i + 1
	cat(i, "/", length(chunks), "\n")
	
	result <- mclapply(chunk, accum.rec, shap, stat, rec, mc.cores=24)
	
	col <- max(unlist(lapply(result, length)))
	mat <- matrix(unlist(result), ncol=col, byrow=TRUE)
	recmap <- rbind(recmap, mat)
	
}

save(recmap, file=sprintf("recmap.random.%s.RData", basename(args[1])))

stop("DONE!")



est <- apply(recmap, 2, mean, na.rm=TRUE)

w <- length(est)
x <- seq(-1.5, 1.5, length.out=w)

plot(x, est, type='l', lwd=1, 
	 xlab="Relative distance to rare variant", ylab="Mean aggregated recombination rate")
x.c <- x[which(x > -0.01 & x < 0.01)]
x.l <- x[which(x < -0.99 & x > -1.01)]
x.r <- x[which(x > 0.99 & x < 1.01)]
polygon(c(x.l, rev(x.l)), c(10000,10000,-1,-1), border=NA, col="grey")
polygon(c(x.r, rev(x.r)), c(10000,10000,-1,-1), border=NA, col="grey")
polygon(c(x.c, rev(x.c)), c(10000,10000,-1,-1), border=NA, col="grey")
lines(x, est, lwd=2)
box(lwd=2)



