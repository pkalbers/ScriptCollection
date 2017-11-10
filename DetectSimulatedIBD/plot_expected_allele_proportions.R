


load("history.RData")

se <- function(x) sqrt(var(x)/length(x))


G = load.bigmatrix(G)
G = G[, ]
G = as.data.frame(G)


frq = rowSums(G) # / (2 * ncol(G))

brk = cut(frq, breaks = seq(0, 2 * ncol(G), length.out = (2 * ncol(G)) + 1), include.lowest = T)

af = 1:(2 * ncol(G))

size = ncol(G)

m0 = matrix(0, nrow = length(af), ncol = size)
m1 = matrix(0, nrow = length(af), ncol = size)
m2 = matrix(0, nrow = length(af), ncol = size)

for (i in 1:size) {
	g = G[, i]

	g = split(g, brk)
	
	m0[, i] = sapply(g, function(g) { if (length(g) == 0) return(0); x = length(which(g == 0)); x / length(g) })
	m1[, i] = sapply(g, function(g) { if (length(g) == 0) return(0); x = length(which(g == 1)); x / length(g) })
	m2[, i] = sapply(g, function(g) { if (length(g) == 0) return(0); x = length(which(g == 2)); x / length(g) })
}

d = rbind(data.table(f = af / (2 * ncol(G)), g = "0", p = apply(m0, 1, mean), e = apply(m0, 1, se)), 
					data.table(f = af / (2 * ncol(G)), g = "1", p = apply(m1, 1, mean), e = apply(m1, 1, se)), 
					data.table(f = af / (2 * ncol(G)), g = "2", p = apply(m2, 1, mean), e = apply(m2, 1, se)))

q = af / max(af)
p = 1 - q

e = rbind(data.table(f = af / (2 * ncol(G)), g = "0", p = p^2), 
					data.table(f = af / (2 * ncol(G)), g = "1", p = 2 * p * q), 
					data.table(f = af / (2 * ncol(G)), g = "2", p = q^2))


ggplot(data = d) + 
	geom_pointrange(aes(x=f, y=p, ymax=p+e, ymin=p-e, colour=g)) + 
	geom_line(data=e, aes(x=f, y=p, linetype=g), colour = "black") +
	theme_bw()




