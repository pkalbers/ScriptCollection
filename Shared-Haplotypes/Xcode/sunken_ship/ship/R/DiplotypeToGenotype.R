
h0_ <- -1:9
h1_ <- -1:9
ph_ <- c("true", "false")


as.binmode <- function(i, n)
	paste(rev(as.numeric(intToBits(i))[1:n]), collapse = '')


M <- c()
H0 <- H1 <- PH <- c()
X <- c()

i <- 0
for (ph in ph_)
	for (h0 in h0_)
		for (h1 in h1_)
		{
			M <- rbind(M, data.frame(ph=ph, h0=h0, h1=h1))
			H0 <- c(H0, sprintf("%s", as.binmode(h0 + 1, 4)))
			H1 <- c(H1, sprintf("%s", as.binmode(h1 + 1, 4)))
			x <- as.character(as.hexmode(i))
			X <- c(X, if (nchar(x) == 1) sprintf("0x0%s", x) else sprintf("0x%s", x))
			i <- i + 1
		}

stop()

M <- data.frame(M, H0=H0, H1=H1, bin=sprintf("0b%s%s", H0, H1), stringsAsFactors = FALSE)


for (i in 1:nrow(M))
{
	cat(M$bin[i], ", ", sep="")
	if (i %% 16 == 0) cat("\n")
}


for (i in 1:nrow(M))
{
	m <- M[i, ]
	cat(sprintf("{%s, %s}, ", 
				if (nchar(as.character(m$h0)) == 1) sprintf(" %d", m$h0) else sprintf("%d", m$h0), 
				if (nchar(as.character(m$h1)) == 1) sprintf(" %d", m$h1) else sprintf("%d", m$h1)
				))
	if (i %% 16 == 0) cat("\n")
}


for (i in 1:length(X))
{
	cat(sprintf("%s, ", X[i]))
	if (i %% 16 == 0) cat("\n")
}



