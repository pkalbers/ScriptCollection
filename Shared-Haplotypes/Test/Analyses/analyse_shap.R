

# args <- commandArgs(TRUE)
# read.table(args[1], TRUE)->stat
# read.table(args[2], TRUE)->shap


shap.files <- dir(pattern="\\.shap$")
stat.files <- dir(pattern="\\.stat$")

if (length(shap.files) != length(stat.files))
	stop("!!!")

shap <- c()
stat <- c()

for(i in 1:length(shap.files)) {
	shap.file <- grep(sprintf("\\.chr%d\\.", i), shap.files, value=TRUE)
	stat.file <- grep(sprintf("\\.chr%d\\.", i), stat.files, value=TRUE)
	
	if (length(shap.file) == 0 || length(stat.file) == 0)
		stop(i)
	
	cat(shap.file, "\n")
	shap.data <- read.table(shap.file, TRUE, stringsAsFactors=FALSE)
	cat(stat.file, "\n")
	stat.data <- read.table(stat.file, TRUE, stringsAsFactors=FALSE)
	
	cat("\n")
	shap <- rbind(shap, shap.data)
	stat <- rbind(stat, stat.data)
}



# physical length
phylength.distr <- function(shap) {	
	splt <- split(shap$length / 1000, shap$mac) # in Kb
	boxplot(splt, outline=FALSE, xlab="Minor allele count (MAC)", ylab="Physical length [Kb]", lwd=2, cex.axis=0.9)
	box(lwd=2)
	
	x <- as.numeric(names(splt))
	y <- unlist(lapply(splt, mean))
	plot(x, y, xlab="Minor allele count (MAC)", ylab="Mean physical length of segments", lwd=3, xaxt="n", type='l')
	axis(side=1, at=as.numeric(names(splt)), labels=names(splt), cex.axis=0.9)
	box(lwd=2)
}

# genetic length
genlength.distr <- function(shap) {
	splt <- split(shap$n_variants / 1000, shap$mac) # in Kb
	boxplot(splt, outline=FALSE, xlab="Minor allele count (MAC)", ylab="Genetic length [Kb]", lwd=2, cex.axis=0.9)
	box(lwd=2)
}

# singletons
singleton.distr <- function(shap) {
	splt <- split(shap$n_singeltons, shap$mac) 

	x <- as.numeric(names(splt))
	y <- unlist(lapply(splt, mean))
	plot(x, y, xlab="Minor allele count (MAC)", ylab="Mean number of singletons per segment", lwd=3, xaxt="n", type='l', ylim=c(0, max(y)))
	axis(side=1, at=as.numeric(names(splt)), labels=names(splt), cex.axis=0.9)
	box(lwd=2)
	

	splt <- split(shap$n_singeltons / shap$n_variants, shap$mac) 
	
	
	x <- as.numeric(names(splt))
	y <- unlist(lapply(splt, mean))
	plot(x, y, xlab="Minor allele count (MAC)", ylab="Mean singleton proportion per segment", 
		 lwd=2, ylim=c(0, max(y)), xaxt="n", type='l')
	axis(side=1, at=x, labels=as.character(x), cex.axis=0.9)
	box(lwd=2)
	
	
	
	
	lambda <- fitdistr(y, "exponential")$estimate
	theta  <- 1 / lambda
	lines(x, c(y[1], y[-1]+ diff(y[1] * exp(-(x / theta)))), col="red", lwd=2)
	
	
	
	
	barplot(y, names.arg=names(splt), border=NA, xpd=FALSE, axes=TRUE, cex.names=0.9,
			xlab="Minor allele count (MAC)", ylab="Number of singletons")
	
	
	splt <- split(shap$n_singeltons, shap$mac)
	
	singl.sum <- lapply(splt, sum)
	singl.frc <- lapply(splt, function(x) sum(x / 6724311)) # = num. singletons
	
	boxplot(splt, outline=FALSE, xlab="Minor allele count (MAC)", ylab="Fraction of singletons per segment", lwd=2, cex.axis=0.9)
	
	box(lwd=2)
}



# split() ...
x <- as.numeric(names(splt))
y <- unlist(lapply(splt, length))
plot(x, y, xlab="Minor allele count (MAC)", ylab="Total number of segments", lwd=3, xaxt="n", type='l')
axis(side=1, at=as.numeric(names(splt)), labels=names(splt), cex.axis=0.9)
box(lwd=2)

shap.num <- unlist(lapply(splt, length))
plot(shap.num, names.arg=names(shap.num), border=NA, xpd=FALSE, axes=TRUE, cex.names=0.9,
		xlab="Minor allele count (MAC)", ylab="Number of pair-wise shared haplotypes")



# Kimura Ota equations for neutral allele age
KO.eq.naa <- function(x)
	( ( -2 * x ) / ( 1 - x ) ) * log(x)

# mutation rate
u <- 1e-08
Ne <- 20000
t <- 4 * Ne * u

x <- sort(unique(stat$maf))
x <- seq(2, 2657, length.out=1000) / (2657 * 2)
age <- KO.eq.naa(x)
plot(x, age, type="l", ylim=c(max(age), min(age)), bty='n', lwd=2, col="grey")

