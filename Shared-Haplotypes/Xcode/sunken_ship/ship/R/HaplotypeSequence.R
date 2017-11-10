#
# make look up tables for compressed genotype data
#

states <- c('A', 'B')

M <- data.frame()
N <- c()

for (i in 1:length(states)) {
	for (j in 1:length(states)) {
		for (k in 1:length(states)) {
			for (l in 1:length(states)) {
				for (m in 1:length(states)) {
					for (n in 1:length(states)) {
						for (o in 1:length(states)) {
							for (p in 1:length(states)) {
				
								comb <- c(i, j, k, l, m, n, o, p) - 1
								code <- sprintf("{%s}", paste(sprintf("%s", states[c(i, j, k, l, m, n, o, p)]), collapse = ", "))
								
								M <- rbind(M, as.list(comb))
								N <- c(N, code)
							}
						}
					}
				}
			}	
		}	
	}
}


print(t(matrix(sprintf("%d,", 0:(length(N)-1)), ncol = 8)), quote = FALSE, row.names = FALSE)

print(sprintf("%s,", N), quote = FALSE, row.names = FALSE)


codes <- sprintf("%s, ", N)

for (i in 1:length(codes)) {
	cat(codes[i])
	if (i %% 2 == 0) cat("\n")
}

hex <- sprintf("0x%s, ", as.hexmode(0:255))
for (i in 1:length(hex)) {
	cat(hex[i])
	if (i %% 16 == 0) cat("\n")
}




