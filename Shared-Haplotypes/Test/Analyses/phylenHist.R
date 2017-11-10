

#shap <- read.table("data.GoT2D.chr20.paper_integrated_snps_indels_sv_beagle_thunder.shap", TRUE, stringsAsFactors=FALSE)

#load("data.GoT2D.chr20.paper_integrated_snps_indels_sv_beagle_thunder.shap.RData")


splt <- split(shap, shap$mac)

phylen <- lapply(splt, function(mac) mac$length)



max <- max(unlist(lapply(phylen, max))) / 1000

png("plot.phylenHist.png", width=8000, height=8000, res=300, pointsize=18)
layout(matrix(1:length(phylen), nrow=5, byrow=TRUE))
mar <- par("mar")

fit <- list()
i <- 0
for (tag in names(phylen)) {
	i <- i + 1
	mac <- as.numeric(tag)
	len <- phylen[[tag]]# / 1000
	med <- median(len / 1000)
	
	fit[[tag]] <- fitdistr(len, "exponential")
	
	par(mar=c(mar[1], mar[2]/2, mar[3]/2, mar[4]))
	hist(len / 1000, freq = FALSE, breaks = seq(0, max, length.out=101), 
		 xlim=c(0, max/3), ylim=c(0, 8e-04), xlab="", ylab="", main="", lwd=1,
		 border="white", col="grey")
	curve(dexp(x, rate = fit[[tag]]$estimate * 1000), col="red", add = TRUE, lwd=2)
	abline(v=med, lwd=2, lty="dashed", col="black")
	
	title(paste("MAC = ", mac), adj=0, line=1, cex=0.5)
	mtext(bquote(paste(lambda, " = ", .(fit[[tag]]$estimate))), adj=1, line=0.5, cex=0.8)
	text(x=round(med, 1), y=8e-04, adj=0, labels=paste("", med))
	
	par(mar=mar)
}

dev.off()


y <- unlist(lapply(fit, function(f) f$estimate))
x <- as.numeric(names(phylen))

mar <- par("mar")
par(mar=c(mar[1], mar[2]*1.5, mar[3], mar[4]))
plot(x, y, xaxt='n', xlab="Minor allele count (MAC)", ylab=bquote(paste("Estimated exponential rate parameter (", lambda, ")")), type='l', lwd=2)
axis(side=1, at=x, labels=as.character(x), cex.axis=0.75)
box(lwd=2)
par(mar=mar)












png("plot.phylenHist.png", width=8000, height=8000, res=300, pointsize=18)
layout(matrix(1:length(phylen), nrow=5, byrow=TRUE))
mar <- par("mar")

for (tag in names(phylen)) {
	mac <- as.numeric(tag)
	pl <- phylen[[tag]]
	
	med <- median(pl)
	
	fit <- fitdistr(pl / chr, "exponential")
	
	par(mar=c(mar[1], mar[2]/2, mar[3]/2, mar[4]))
	hist(pl / chr, freq = FALSE, breaks = seq(0, 1, length.out=101), xlim=c(0, 0.5), ylim=c(0, 70), xlab="", ylab="", main="", lwd=1)
	curve(dexp(x, rate = fit$estimate), col = "#FF000080", add = TRUE, lwd=2)
	abline(v=med / chr, lwd=2, lty="dashed", col="grey")
	
	title(paste("MAC = ", mac), adj=0, line=1, cex=0.5)
	mtext(paste("rate = ", round(fit$estimate, 2)), adj=1, line=0.5, cex=0.8)
	text(x=med / chr, y=70, adj=0, labels=paste("", med))
	
	par(mar=mar)
}

dev.off()


