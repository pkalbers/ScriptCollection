
h0_ <- -1:6
h1_ <- -1:6
ph_ <- 0:1


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
			M <- rbind(M, c(ph=ph, h0=h0, h1=h1))
			H0 <- c(H0, sprintf("%s", as.binmode(h0 + 1, 3)))
			H1 <- c(H1, sprintf("%s", as.binmode(h1 + 1, 3)))
			PH <- c(PH, sprintf("%s", as.binmode(ph, 1)))
			X <- c(X, sprintf("0x%s", as.hexmode(i)))
			i <- i + 1
		}


M <- data.frame(M, PH=PH, H0=H0, H1=H1, bin=sprintf("0b0%s%s%s", PH, H0, H1), stringsAsFactors = FALSE)


for (i in 1:nrow(M))
{
	cat(M$bin[i], ", ", sep="")
	if (i %% 8 == 0) cat("\n")
}


for (i in 1:nrow(M))
{
	m <- M[i, ]
	cat(sprintf("{%s, %s, %s}, ", 
				if (m$h0 == -1) "-1" else sprintf(" %d", m$h0), 
				if (m$h1 == -1) "-1" else sprintf(" %d", m$h1), 
				if (m$ph == 1) " true" else "false"))
	if (i %% 8 == 0) cat("\n")
}




