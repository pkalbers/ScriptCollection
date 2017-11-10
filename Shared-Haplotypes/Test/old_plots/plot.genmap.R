

m <- read.table("summary.marker", header = TRUE, stringsAsFactors = FALSE)


m <- m[ -(which(m$position < 100000)) , ]

m <- m[ -(which(m$position > m$position[1] + 500000)) , ]

write.table(data.frame(pos=m$position, map_rate=m$map_rate, map_dist=m$map_dist, map_source=m$map_source), file = "mapinterp.txt", quote = FALSE, row.names = FALSE)


plot(m$position, m$map_dist, pch='.', xlab="Physical position [1Mb region]", )
points(m$position[which(m$map_source == 'm')], m$map_dist[which(m$map_source == 'm')], col='red', pch='|')
