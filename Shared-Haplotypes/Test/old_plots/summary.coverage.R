

library("parallel")


cov <- read.table("summary.coverage.txt", header = TRUE, stringsAsFactors = FALSE)
cov <- cov[order(paste(cov$sample_group, cov$sample_pop, cov$sample_id, sep="_")),]

cov$sample_id.1 <- NULL
cov$sample_key <- NULL
cov$struct <- NULL
cov$from_pos <- NULL
cov$to_pos <- NULL

splt <- split(cov, cov$sample_id)

x <- max(cov$to_id) + 1


f.select <- function(data, f, cum)
{
	if (cum) {
		
		return(which(data$n_sample <= f))
		
	} else {
		
		return(which(data$n_sample == f))
		
	}
}


res <- mclapply(splt, function(s, x) { 

	n <- split(s, s$n_sample)
	
	fs <- 2:25
	
	res <- data.frame(f=fs, percent.cum=NA, percent.exact=NA) #min.cum=NA, min.exact=NA, max.cum=NA, max.exact=NA, mean.cum=NA, mean.exact=NA, median.cum=NA, median.exact=NA)
	
	for (f in fs)
	{
		fi <- which(res$f == f)
		
		
		sub <- s[f.select(s, f, TRUE), ]
		
		v <- rep(0, x)
		
		for (i in 1: nrow(sub))
		{
			from <- sub$from_id[i] + 1
			to   <- sub$to_id[i] + 1
			
			v[from:to] <- v[from:to] + 1
		}
		
	#	res$min.cum[fi] = min(v)
	#	res$max.cum[fi] = max(v)
	#	res$mean.cum[fi] = mean(v)
	#	res$median.cum[fi] = median(v)
	
	res$percent.cum[fi] = (1 - (length(which(v == 0)) / length(v))) * 100
		
		
		sub <- s[f.select(s, f, FALSE), ]
		
		v <- rep(0, x)
		
		for (i in 1: nrow(sub))
		{
			from <- sub$from_id[i] + 1
			to   <- sub$to_id[i] + 1
			
			v[from:to] <- v[from:to] + 1
		}
		
	#	res$min.exact[fi] = min(v)
	#	res$max.exact[fi] = max(v)
	#	res$mean.exact[fi] = mean(v)
	#	res$median.exact[fi] = median(v)
	
	res$percent.exact[fi] = (1 - (length(which(v == 0)) / length(v))) * 100
	}
	
	return(res)
	
}, x, mc.cores = 25)



############



m <- read.table("table.coverage.cum.txt", header = TRUE, stringsAsFactors = FALSE)
m <- m[order(m$id), ] 
m <- m[order(m$pop), ]
m <- m[order(m$grp), ]


splt <- split(m, m$grp)

res <- list()

for (tag in names(splt))
{
	fs <- paste("f", 2:25, sep="")
	
	x <- c()
	
	for (f in fs)
	{
		x <- c(x, median(splt[[tag]][[f]]))
	}
	
	names(x) <- fs
	
	res[[tag]] <- x
}

data <- c()
for (tag in names(res))
	data <- rbind(data, res[[tag]])
data <- as.data.frame(data)
data <- cbind(grp=names(res), data)

write.table(data, file = "table.combined.cum.grp.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)



res <- list()

for (tag in names(splt))
{
	fs <- paste("f", 2:25, sep="")
	
	sub <- split(splt[[tag]], splt[[tag]]$pop)
	
	for (s in names(sub))
	{
		x <- c()
		
		for (f in fs)
		{
			x <- c(x, median(sub[[s]][[f]]))
		}
		
		names(x) <- fs
		
		res[[s]] <- x
	}
}

data <- c()
for (tag in names(res))
	data <- rbind(data, res[[tag]])
data <- as.data.frame(data)
data <- cbind(grp=names(res), data)

write.table(data, file = "table.combined.cum.pop.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)




w <- m[which(m$f25 == 100), ]

w <- cbind(tag=paste(w$grp, w$pop, w$id, sep=" "), w)

w$grp <- NULL
w$pop <- NULL
w$id <- NULL

write.table(w, file = "table.combined.cum.100.txt", quote = TRUE, row.names = FALSE, col.names = TRUE)




###########



d <- read.table("table.coverage.cum.txt", header = TRUE, stringsAsFactors = FALSE)

pop <- d$pop
grp <- d$grp
id <- d$id

d$grp <- NULL
d$pop <- NULL
d$id <- NULL

col <- c()
color <- rainbow(length(unique(grp)))
names(color) <- unique(grp)

for (i in 1:length(grp))
	col <- c(col, color[which(names(color) == grp[i])])

m<-(as.matrix(d))

rownames(m) <- id


pca <- prcomp(m, retx=TRUE, center=TRUE, scale.=TRUE) 

percent <- round((((pca$sdev)^2 / sum(pca$sdev^2))*100)[1:2])
loadings <- pca$rotation
rownames(loadings) <- colnames(m)
scores <- pca$x

PCA1 <- scores[,1]
PCA2 <- scores[,2]


Group <- grp

q <- qplot(PCA2, PCA1, main="PCA of chromosomeal coverage, f cumulative", geom="blank", xlab = paste("PCA2 (", percent[2], "%)", sep = ""), ylab = paste("PCA1 (", percent[1], "%)", sep = ""))
q <- q + geom_point(aes(colour = Group), size = 3, alpha = 3/4)

ggsave(filename = "table.coverage.cum.PCA.png", width = 10, height = 10)

