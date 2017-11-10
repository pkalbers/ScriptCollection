
args <- commandArgs(TRUE)

read.table(args[1], TRUE)->shap


library("parallel")



nested.shap <- function(i, shap.from, shap.to, shap.mac, shap.i0, shap.i1) {
	from <- shap.from[i]
	to   <- shap.to[i]
	mac  <- shap.mac[i]
	
	nested <- which(shap.mac < mac & shap.from >= from & shap.to <= to)
	
	if (length(nested) == 0)
		return(NULL)
	
	i0   <- shap.i0[i]
	i1   <- shap.i1[i]
	
	s0 <- shap.i0[nested]
	s1 <- shap.i1[nested]
	
	sam.a <- which(s0 == i0 & s1 == i1)
	sam.b <- which(s0 == i1 & s1 == i0)
	sam   <- c(sam.a, sam.b)
	
	if (length(sam) == 0)
		return(NULL)
	
	macs <- shap.mac[nested][sam]
	macs <- tabulate(macs, max(unique(shap.mac)))
	
	return(c(x=mac, macs))
}



full <- 1:nrow(shap)
until <- 50000
chunks <- split(full, ceiling(seq_along(full)/until))

nest <- c()

i <- 0
for (chunk in chunks) {	
	i <- i + 1
	cat(i, "/", length(chunks), "\n")
	
	result <- mclapply(chunk, nested.shap, shap$from, shap$to, shap$mac, shap$i0, shap$i1, mc.cores=3)
	
	width <- max(unlist(lapply(result, length)))
	
	matr <- matrix(unlist(result), ncol=width, byrow=TRUE)
	nest <- rbind(nest, matr)
	
	save(nest, file=sprintf("nest.%s.RData", basename(args[1])))
}

stop("DONE")


nest <- as.data.frame(nest)
macs <- split(nest[, -1], nest[, 1])

est <- lapply(macs, function(m) apply(m, 2, sum))
err <- lapply(macs, function(m) apply(m, 2, function(x) (sd(x) / sqrt(length(x)))))
for (name in names(est)) {
	max <- as.numeric(name) - 1
	est[[name]] <- est[[name]][2:max]
	err[[name]] <- err[[name]][2:max]
	names(est[[name]]) <- 2:max
	names(err[[name]]) <- 2:max
}

max <- max(unlist(lapply(est, max)))


png(filename="nested.chr20.png", height=30000, width=10000)
h <- (1:length(est))+2
h[1] <- h[1]*0.75
layout(matrix(1:length(est), length(est), 1, byrow = TRUE), heights=h)
for (name in names(est)) {
	mai <- par("mai")
	par(mai=c(mai[1], mai[2]+1, mai[3:4]))
	y <- barplot(rev(est[[name]]), border=NA, xlim=c(0, max), xaxt='n', yaxt='n', horiz=TRUE, xpd=FALSE)
	labs <- rev(names(est[[name]]))
	text(cex=10, y=y, x=rep(-0.02, length(y)), labs, xpd=TRUE)
	par(mai=mai)
}
dev.off()




cols <- sort(as.numeric(unique(unlist(lapply(est, names)))))
rows <- as.numeric(names(est))

tab <- matrix(NA, length(rows), length(cols))
colnames(tab) <- as.character(cols)
rownames(tab) <- as.character(rows)
for (row in rows) {
	dat <- est[[as.character(row)]]
	for (col in cols) {
		if (as.character(col) %in% names(dat))
			tab[as.character(row), as.character(col)] <- dat[as.character(col)]
	}
}

max <- sum(tab, na.rm=TRUE)
tab <- (tab / max) * 100


