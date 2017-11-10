


library(data.table)


for (i in 1:22) {
	cat(i, "\n")
	
	load(sprintf("pruned.chr%d.net.RData", i))
	
	sub = net[sample(1:nrow(net), size = 10000, replace = F), ] 
	
	cat(sample(sub$Position), sep = "\n", file = sprintf("selected_chr%d.neutral.txt", i))
	
}




