

args <- commandArgs(TRUE)
shap.file <- args[1]

shap <- read.table(shap.file, TRUE, stringsAsFactors=FALSE)




cantor <- function(i, j) 
	0.5 * (i+j) * (i+j+1) + j



B <- U <- L <- matrix(0, 2657, 2657)
str <- paste(shap$i0 + 1, j=shap$i1 + 1)
num <- table(str)
ind <- strsplit(names(num), " ", fixed=TRUE)
i0 <- unlist(lapply(ind, function(x) as.numeric(x[1])))
i1 <- unlist(lapply(ind, function(x) as.numeric(x[2])))

for (x in 1:length(num)) {
	i <- i0[x]
	j <- i1[x]
	n <- as.numeric(num[x])
	B[i, j] <- n
	B[j, i] <- n
	U[i, j] <- n
	L[j, i] <- n
}

save(B, U, L, file="varsharing.allChr.RData")

x <- B#[1:1000, 1:1000]
#x <- (1 - (x / max(x)))^10
x <- log10(x)
x[which(is.infinite(x))] <- 0

png("popstruct.log10.allChr.png", width=4000 * 1.5, height=3800 * 1.5, res=300, pointsize=18)

q <- qplot(x=Var1, y=Var2, data=melt(x), fill=value, geom="tile", xlab="", ylab="")
q+theme_bw()+scale_fill_gradient(low="white", high="black")

dev.off()




x <- B[1500:2000, 1500:2000]
x <- exp((x / max(x)))
png("popstruct.exp.allChr.png", width=4000, height=3900, res=300, pointsize=18)
q <- qplot(x=Var1, y=Var2, data=melt(x), fill=value, geom="tile")
q+theme_bw()+scale_fill_gradient(low="white", high="black")
dev.off()




heatmap((x / max(x)), symm=TRUE, Rowv=NA)



P <- M / max(M)



index <- cantor(shap$i0 + 1, shap$i1 + 1)

count <- table(index)

df <- data.frame(index=as.numeric(names(count)), count=count)

fill.M <- function (i, index) {
	row <- rep(0, 2657)
	for (j in 1:2657) {
		want <- cantor(i, j)
		find <- length(which(index == want))
		row[j] <- find
	}
	return(row)
}


result <- mclapply(1:2657, fill.M, index, mc.cores=20)
#result <- mclapply(1:3, fill.M, index, mc.cores=3)

M <- Reduce(rbind, Result)
row.names(M) <- NULL

save(M, file=sprintf("varsharing.%s.RData", shap.file))


