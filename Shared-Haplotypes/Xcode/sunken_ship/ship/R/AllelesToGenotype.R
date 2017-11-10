

a <- 0:21

L<-list()

for (i in a) 
	for (j in a) 
		L<-c(L, list(sort(c(i, j))))

L <- unique(L)


C <- unlist(lapply(L, function(x) min(x) + choose(max(x)+1, 2)))

for (i in 1:length(C)) {
	L[[i]] <- c(L[[i]], C[i])
}




M <- matrix(NA, length(a), length(a))
for (l in L) M[l[1]+1, l[2]+1] <- l[3]
for (l in L) M[l[2]+1, l[1]+1] <- l[3]
apply(M, 1, function(x) cat(paste(sprintf("0x%s", as.hexmode(x)), collapse=", "), ",\n", sep=""))


N <- matrix(NA, length(L), 2)
for (l in L)  { N[l[3]+1, 1] <- l[1]; N[l[3]+1, 2] <- l[2]; }
for (i in 1:nrow(N)) {
	if (N[i, 1] == 0) cat("\n")
	cat(sprintf("{%d, %d}, ", N[i, 1], N[i, 2]))
}

