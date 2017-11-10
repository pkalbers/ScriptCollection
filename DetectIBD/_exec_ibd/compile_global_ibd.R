#
# time & freq dependend genotype pair emissions
#

library(data.table)


args = commandArgs(T)

prefix = "OutOfAfricaHapMap20"
ibd.file = args[1]

cat("Loading...\n")
load(sprintf("../data.%s.RData", prefix))
load(sprintf("../data.%s.G.RData", prefix)) # genotypes
load(ibd.file)
cat("...OK\n")


cut.time = function(t, brk = c(0, 10^(0:6)), lab = sprintf("% 8d", c(10^(0:6)))) {
	as.character(cut(t, breaks = brk, labels = lab, include.lowest = T))
}

cut.geno = function(g0, g1, mask = c("00"="00", "01"="01", "10"="01", "02"="02", "20"="02", "11"="11", "12"="12", "21"="12", "22"="22")) {
	g = sprintf("%d%d", g0, g1)
	as.vector(mask[g])
}

make.matrix = function(len = nrow(G)) {
	m = matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
	M = list()
	for (tag in sprintf("% 8d", c(10^(0:6)))) {
		M[[tag]] = m
	}
	M
}


count.pairs = function(ibd, pair, min.size = 10, mask = c("00"=1, "01"=2, "02"=3, "11"=4, "12"=5, "22"=6)) {
	count = make.matrix()
	
	ibd$n = ibd$rhs.index - ibd$lhs.index
	
	del = which(ibd$n < min.size)
	if (length(del) > 0) {
		ibd = ibd[-del, ]
	}
	
	time = cut.time(ibd$time)
	
	for (i in 1:nrow(ibd)) {
		t = time[i]
		
		lhs = ibd$lhs.index[i] + 2
		rhs = ibd$rhs.index[i] - 2
		
		rng = lhs : rhs
		
		coord = matrix(c(rng, mask[pair[rng]]), ncol = 2, byrow = F)
		count[[ t ]][ coord ] = count[[ t ]][ coord ] + 1
	}
	
	count
}


###


cat("Observed genotypes per site ...\n")

count = make.matrix()

i = 0
for (pair in names(ibd)) {
	i = i + 1
	cat(sprintf("%d: %s\n", i, pair))
	
	cmb = ibd[[pair]]
	
	g0 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\1", pair))
	g1 = as.numeric(sub("^([0-9]+) ([0-9]+)$", "\\2", pair))
	g0g1 = cut.geno(G[, g0], G[, g1])
	
	for (tag in names(cmb)) {
		tmp = count.pairs(cmb[[tag]], g0g1)
		for (time in names(tmp)) {
			count[[time]] = count[[time]] + tmp[[time]]
		}
	}
}
cat("\n")


save(count, file = sprintf("result.%s", ibd.file))





###
stop("DONE")
###


library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20"

load(sprintf("../data.%s.RData", prefix))



cut.freq = function(f, brk = (0:101)/101, lab = sprintf("%.3f", (1:101)/101)) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

make.matrix = function(len = length(POS)) {
	m = matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
	M = list()
	for (tag in sprintf("% 8d", c(10^(0:6)))) {
		M[[tag]] = m
	}
	M
}


res.files = dir(pattern = "^result\\..+\\.RData$")

tmp = make.matrix()
for (res.file in res.files) {
	cat(res.file, "\n")
	
	load(res.file)
	
	for (tag in names(count)) {
		tmp[[tag]] = tmp[[tag]] + count[[tag]]
	}
}
count = tmp
rm(tmp)


count$`      10` = count$`      10` + count$`       1`
count$`       1` = NULL

count$`  100000` = count$`  100000` + count$` 1000000`
count$` 1000000` = NULL


for (i in 2:length(count)) {
	count[[i]] = count[[i]] + count[[i-1]]
}


cat("Observed genotypes per frequency bin ...\n")

freq = list()

for (time in names(count)) {
	tmp = as.data.frame(count[[time]])
	tmp = split(tmp, cut.freq(AAF))
	tmp = lapply(tmp, function(x) {
		x = colSums(x)
		x = x / sum(x)
		as.data.table(as.list(x))
	})
	tag = names(tmp)
	tmp = rbindlist(tmp)
	rownames(tmp) = tag
	
	freq[[time]] = tmp
}
cat("\n")



### plotting ...


col = c("00" = "#6082E5",
				"01" = "#46B29D",
				"02" = "#B571E3",
				"11" = "#DBAC12",
				"12" = "#969130",
				"22" = "#DE473A")


d = NULL

for (t in names(freq)) {
	f = as.numeric(rownames(freq[[t]]))
	
	for (p in colnames(freq[[t]])) {
		x = freq[[t]][[p]]
		
		d = rbind(d, data.table(time = t, 
														freq = f, 
														pair = p,
														prop = x))
	}
}

del = which(is.na(d$prop))
if (length(del) > 0) {
	d = d[-del, ]
}


q = ((0:100)/100)
p = 1 - q

expect.non = rbind(data.table(freq = q, pair = "00", prop = p^4),
									 data.table(freq = q, pair = "01", prop = 4*(p^3)*q),
									 data.table(freq = q, pair = "02", prop = 2*(p^2)*(q^2)),
									 data.table(freq = q, pair = "11", prop = 4*(p^2)*(q^2)),
									 data.table(freq = q, pair = "12", prop = 4*p*(q^3)),
									 data.table(freq = q, pair = "22", prop = q^4))
expect.non$type = "Expected under NON"

expect.ibd = rbind(data.table(freq = q, pair = "00", prop = p^3),
									 data.table(freq = q, pair = "01", prop = 2*(p^2)*q),
									 data.table(freq = q, pair = "02", prop = 0),
									 data.table(freq = q, pair = "11", prop = ((p^2)*q) + (p*(q^2))),
									 data.table(freq = q, pair = "12", prop = 2*p*(q^2)),
									 data.table(freq = q, pair = "22", prop = q^3))
expect.ibd$type = "Expected under IBD"

e = rbind(expect.non, expect.ibd)
r = cbind(expect.non, ibd = expect.ibd$prop)




d$time[which(d$time == "      10")] = "    t ≤ 10    "
d$time[which(d$time == "     100")] = "   t ≤ 100   "
d$time[which(d$time == "    1000")] = "  t ≤ 1000  "
d$time[which(d$time == "   10000")] = " t ≤ 10000 "
d$time[which(d$time == "  100000")] = "t > 10000"



gg = ggplot(d) +
	facet_wrap(~time, nrow = 1) +
	geom_line(data = expect.ibd, aes(x = freq, y = prop, group = pair), linetype = "33", alpha = 2/3) +
	geom_line(data = expect.non, aes(x = freq, y = prop, group = pair), linetype = "11", alpha = 2/3) +
	geom_ribbon(data = r, aes(x=freq, ymin = prop, ymax = ibd, fill = pair), alpha = 1/3) +
	geom_point(aes(x = freq, y = prop, colour = pair), alpha = 0.8, size = 1) +
	coord_cartesian(xlim = c(-0.025, 1.025), ylim = c(-0.01, 1.01), expand = F) +
	scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_colour_manual(values = col) +
	scale_fill_manual(values = col) +
	theme_few() +
	theme(aspect.ratio = 2.8125,
				panel.border = element_rect(fill = NA, colour = "black", size = 1/2),
				legend.title = element_blank(),
				legend.background = element_rect(fill = "white", colour = "black", size = 1/3),
				legend.position = c(1 - (1/10.5),0.825), 
				panel.grid = element_line(colour = "grey90", size = 1/2),
				panel.grid.major = element_line(colour = "grey90", size = 1/2),
				panel.grid.minor = element_blank()) +
	xlab("Allele frequency") + ylab("Genotype pair proportion") +
	guides(color = guide_legend(override.aes = list(alpha = 1)))
gg

ggsave(filename = "./_plot.genpair_tmrca.pdf", plot = gg, width = 11, height = 6.5)


ggplot(d) +
	facet_wrap(~pair) +
	geom_point(aes(x = freq, y = prop, colour = time), size = 0.25) 













