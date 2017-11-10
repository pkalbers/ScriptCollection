

#shap <- read.table("data.GoT2D.chr20.paper_integrated_snps_indels_sv_beagle_thunder.shap", TRUE, stringsAsFactors=FALSE)

load("data.GoT2D.chr20.paper_integrated_snps_indels_sv_beagle_thunder.shap.RData")

library("parallel")

splt <- split(shap, shap$mac)

# del <- which(as.numeric(names(splt)) <= 2)
# if (length(del) != 0)
# 	splt[del] <- NULL


map.phylen.by.mac <- function(mac.group) {
	subs <- split(mac.group, mac.group$position)
	rm(mac.group)
	
	return(unlist(lapply(subs, function(sub) sub$length)))
}

phylen <- mclapply(splt, map.phylen.by.mac, mc.cores=4)


chr <- max(shap$to) - min(shap$from)

png("plot.phylenHist.chr20.png", width=8000, height=8000, res=300, pointsize=18)
layout(matrix(1:length(phylen), nrow=5, byrow=TRUE))
mar <- par("mar")

for (tag in names(phylen)) {
	mac <- as.numeric(tag)
	pl <- phylen[[tag]]
	
	med <- median(pl)
	
	fit <- fitdistr(pl / chr, "exponential")
	
	par(mar=c(mar[1], mar[2]/2, mar[3]/2, mar[4]))
	hist(pl / chr, freq = FALSE, breaks = 50, xlim=c(0, 0.5), ylim=c(0, 70), xlab="", ylab="", main="", lwd=1)
	curve(dexp(x, rate = fit$estimate), col = "#FF000080", add = TRUE, lwd=2)
	abline(v=med / chr, lwd=2, lty="dashed", col="grey")
	
	title(paste("MAC = ", mac), adj=0, line=1, cex=0.5)
	mtext(paste("rate = ", round(fit$estimate, 2)), adj=1, line=0.5, cex=0.8)
	text(x=med / chr, y=70, adj=0, labels=paste("", med))
	
	par(mar=mar)
}

dev.off()







stop("DONE")






lapply(phylen, function(pl, chr) fitdistr(pl/chr, "exponential")$estimate, chr)





map.brkpts.by.mac <- function(mac.group) {
	subs <- split(mac.group, mac.group$position)
	rm(mac.group)
	
	left  <- c()
	right <- c()
	for (sub in subs) {
		pos <- sub$pos[1]
		
		dist.l <- abs(sub$from - pos)
		dist.r <- abs(sub$to   - pos)
		
		#dist.l <- dist.l / max(dist.l)
		#dist.r <- dist.r / max(dist.r)
		
		left  <- c(left,  dist.l * -1)
		right <- c(right, dist.r)
	}
	
	del <- which(left == 0)
	if (length(del) != 0)
		left <- left[-del]
	
	del <- which(right == 0)
	if (length(del) != 0)
		right <- right[-del]
	
	return(list(n=length(subs), left=left, right=right))
}

brkpts <- mclapply(splt, map.brkpts.by.mac, mc.cores=4)



dens <- list()
for (tag in names(phylen)) {
	#dens[[tag]] <- density(c(brkpts[[tag]]$left, brkpts[[tag]]$right))
	dens[[tag]] <- density(phylen[[tag]], kernel="epanechnikov")
}

xlim <- c(min(unlist(lapply(dens, function(d) min(d$x)))), max(unlist(lapply(dens, function(d) max(d$x)))))
ylim <- c(min(unlist(lapply(dens, function(d) min(d$y)))), max(unlist(lapply(dens, function(d) max(d$y)))))

plot(NA, xlim=xlim, ylim=ylim)

i <- 0
for (tag in names(dens)) {
	i <- i + 1
	x <- 
	lines(dens[[tag]], col=terrain.colors(length(dens), alpha=0.9)[i], lwd=2)
}
legend("topright", legend=names(dens), col=terrain.colors(length(dens), alpha=0.9), cex=0.5, lty="solid", lwd=2)
abline(v=0)



plot(as.numeric(names(brkpts)), unlist(lapply(brkpts, function(x) x$n)))


stop("DONE")

splt <- split(shap, shap$position)

library("parallel")

find.internal.break <- function(x, left=TRUE) {
	if (nrow(x) < 2)
		return(NA)
	pos <- x$pos[1]
	dst <- 0
	if (left)
		dst <- abs(x$from - pos)
	else
		dst <- abs(x$to - pos)
	dst <- dst / max(dst)
	return(dst)
}

left  <- unlist(mclapply(splt, find.internal.break, left = TRUE, mc.cores=8))
right <- unlist(mclapply(splt, find.internal.break, left = FALSE, mc.cores=8))


del1 <- which(right == 1)
del2 <- which(is.na(right))
r <- right[-c(del1, del2)]

del1 <- which(left == 1)
del2 <- which(is.na(left))
l <- left[-c(del1, del2)]

dist <- c(l * -1, r)

plot(density(dist, na.rm=TRUE), lwd=2, main="", xlab="Relative distance to rare variant", ylab="Density of breakpoints")
abline(v=0, lwd=5, col="white")
abline(v=-1.01, lwd=8, col="white")
abline(v=1.01, lwd=8, col="white")
abline(h=0, lwd=1, lty="dashed")
abline(v=0, lwd=3, col="grey")
abline(v=-1, lwd=3, col="grey")
abline(v=1, lwd=3, col="grey")
box(lwd=2)


save(left, right, dist, file="plot.breakdist.RData")





