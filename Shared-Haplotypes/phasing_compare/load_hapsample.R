
data.file = "data.sim.5000.RData"

load(data.file)




samfile = "phased.sim.5000.sample"

sam = read.table(samfile, header = FALSE, stringsAsFactors = F)

sam = sam[-(1:2), ]

sam = sam[, 1]

sam = I[sam]

sam = unlist(strsplit(sam, " ", T))

sam = as.vector(sam)


hapfile = "phased.sim.5000.haps"

lines = readLines(hapfile, n = -1)

lines = strsplit(lines, " ", T)

PH = matrix(NA, nrow = length(lines), ncol = 5000, dimnames = list(NULL, sam))

for (i in 1:length(lines)) {
	line = lines[[i]]
	k = as.numeric(sub("^SNP_([0-9]+)$", "\\1", line[2]))
	cat(k, "\n")
	
	line = line[-(1:5)]
	line = as.numeric(line)
	
	PH[i, ] = line
}


save(PH, file = "phased.sim.5000.RData")


