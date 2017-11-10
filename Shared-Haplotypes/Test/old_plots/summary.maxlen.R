

data <- read.table("summary.maxlen.txt", header = TRUE, stringsAsFactors = FALSE)

splt <- split(data, data$n_samples)

num <- unlist(lapply(splt, function(x) length(x$gen_dist)))
phy <- unlist(lapply(splt, function(x) median(x$phy_dist)))
gen <- unlist(lapply(splt, function(x) median(x$gen_dist)))

num <- num[-1]
phy <- phy[-1]
gen <- gen[-1]


write.table(data.frame(f = names(num), num = num, phy = phy, gen = gen), file = "result_summary.maxlen.txt", quote = FALSE, row.names = FALSE)


####


data <- read.table("./biallelic/summary.nodes.txt", header = TRUE, stringsAsFactors = FALSE)

data <- cbind(data, all_nodes = data$l_nodes + data$r_nodes)

splt <- split(data, data$n_samples)

b<-unlist(lapply(splt, function(x) mean(x$all_nodes)))

