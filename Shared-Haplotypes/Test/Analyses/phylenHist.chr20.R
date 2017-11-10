

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


