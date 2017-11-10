
library(data.table)
library(ggplot2)
library(ggthemes)


D = fread("vanilla.marker.txt", header = T)
M = read.table("vanilla.mutations.txt", header=T)

rec = 1e-08

act = D$Position
tru = M$position


del = which(diff(tru) < rec)
if (length(del) > 0) {
	act = act[-del]
	tru = tru[-del]
}


dist = (tru / 1e6) * rec * 100 * 1e6
dist = dist - dist[1]

rate = c(diff(dist) / (diff(act) / 1e6), 0)



map = data.table(Chromosome = "chr1", "Position(bp)" = act, "Rate(cM/Mb)" = sprintf("%.10f", rate), "Map(cM)" = sprintf("%.10f", dist))
			
write.table(map, file = "vanilla.map", append = F, quote = F, sep = "	", row.names = F, col.names = T)


