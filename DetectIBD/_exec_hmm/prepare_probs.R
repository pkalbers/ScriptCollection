#
# Preparing HMM probs
#

library(data.table)
library(ggplot2)
library(ggthemes)


prefix = "OutOfAfricaHapMap20.GenErr_1000G"
prefix = "OutOfAfricaHapMap20"

original = "OutOfAfricaHapMap20"

load(sprintf("data.%s.RData", prefix))

load(sprintf("../data.%s.G.RData", original))
sim = G
rm(G)

load(sprintf("data.%s.G.RData", prefix))



cut.freq = function(f, brk = (0:200)/200, lab = sprintf("%.8f", ((1:200)/200) - (0.5/200))) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}
cut.freq = function(f, brk = (0:500)/500, lab = sprintf("%.8f", ((1:500)/500) - (0.5/500))) {
	cut(f, breaks = brk, labels = lab, include.lowest = T)
}

make.matrix = function(len = length(POS)) {
	matrix(0, nrow = len, ncol = 6, dimnames = list(NULL, c("00", "01", "02", "11", "12", "22")))
}


res.files = dir(pattern = sprintf("^result\\.local_ibd_counts\\.%s\\..+\\.RData$", prefix), path = "./local_ibd/", full.names = T)[1:40]

num = c(IBD=0, NON=0)

ibd = make.matrix()
non = make.matrix()

mibd = NULL
mnon = NULL

for (res.file in res.files) {
	cat(res.file, "\n")
	
	load(res.file)
	
	num["IBD"] = num["IBD"] + nrow(count$IBD)
	num["NON"] = num["NON"] + nrow(count$NON)
	
	ibd = ibd + count$IBD
	non = non + count$NON
	
	mibd = rbind(mibd, count$IBD)
	mnon = rbind(mnon, count$NON)
}
num


# grid plot

mibd = ibd #/ rowSums(ibd)
mnon = non #/ rowSums(non)

mibd = as.data.table(mibd)
mnon = as.data.table(mnon)

mibd$frq = AAF #rep(AAF, 40)
mnon$frq = AAF #rep(AAF, 40)

mibd = split(mibd, mibd$frq)
mibd = lapply(mibd, function(x) { if(nrow(x)==0) return(NULL); z=colSums(x)[1:6]; z = data.table(t(z/sum(z)), frq = x$frq[1]);  z  })
mibd = rbindlist(mibd)

mnon = split(mnon, mnon$frq)
mnon = lapply(mnon, function(x) { if(nrow(x)==0) return(NULL); z=colSums(x)[1:6]; z = data.table(t(z/sum(z)), frq = x$frq[1]);  z  })
mnon = rbindlist(mnon)

mibd = melt.data.table(mibd, measure.vars = colnames(ibd))
mnon = melt.data.table(mnon, measure.vars = colnames(non))

mibd$state = "ibd"
mnon$state = "non"

d = rbind(mnon, mibd)

del = which(is.na(d$value))
if (length(del) > 0) {
	d = d[-del, ]
}

get.delta = function(d, x = 5e6) {
	if (length(unique(d$variable)) != 1) stop()
	if (length(unique(d$state)) != 1) stop()
	gt = d$variable[1]
	if (nrow(d) > x)
		d = d[sample(1:nrow(d), x),]
	tmp = NULL
	q = d$frq
	p = 1- q
	if (d$state[1] == "non") {
		if (gt == "00") tmp = p^4
		if (gt == "01") tmp = 4 * p^3 * q
		if (gt == "02") tmp = 2 * p^2 * q^2
		if (gt == "11") tmp = 4 * p^2 * q^2
		if (gt == "12") tmp = 4 * p * q^3
		if (gt == "22") tmp = q^4
	}
	if (d$state[1] == "ibd") {
		if (gt == "00") tmp = p^3
		if (gt == "01") tmp = 2 * p^2 * q
		if (gt == "02") tmp = 0
		if (gt == "11") tmp = (p^2 * q) + (p * q^2)
		if (gt == "12") tmp = 2 * p * q^2
		if (gt == "22") tmp = q^3
	}
	d$expected = tmp
	d$delta = d$value - tmp
	d
}

# d$cat = round(d$frq * 1000)
# p = split(d, list(d$cat, d$variable, d$state))


p = split(d, list(d$variable, d$state))
p = rbindlist(lapply(p, get.delta))

p$state = factor(p$state, levels = c("non", "ibd"), ordered = T)

p$cat = round(p$frq * 100) / 100

q = p[which(p$frq == 0 | p$frq == 1),]
for (i in 1:1000) { p = rbind(p, q) }

gg = ggplot(p) +
	facet_grid(state~variable) +
	#geom_rect(data = data.table(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey80") +
	stat_bin2d(aes(frq, delta), binwidth = c(0.005, 0.001), geom = "raster") +
	#stat_summary(aes(cat, delta), fun.y = "mean", geom = "line", color = "black", alpha = 0.75) +
	stat_summary(aes(cat, delta), fun.y = "mean", geom = "line", color = "black", size = 1/2) +
	geom_hline(yintercept = 0, color = "black", size=1/2, linetype="11") +
	#geom_density_2d(aes(frq, delta), bins = 5, size = 1/3, color = "black", alpha = 0.8) +
	# scale_fill_gradientn(colours = c("white", "yellow", "orange", "orangered", "red", "darkred"), na.value = "black", limits = c(1, 100), trans = "log10",
	# 										 breaks = (10^(0:8)), labels = sprintf("%d", 10^(0:8))) +
	scale_fill_gradientn(colours = c("grey65", "grey60"), na.value = "black", #limits = c(0.0001, 1), trans = "log10",
											 breaks = (10^(0:8)), labels = sprintf("%d", 10^(0:8))) +
	scale_x_continuous(breaks = (0:4)/4, labels = c(" 0", "0.25", "0.5", "0.75", "1 ")) +
	scale_y_continuous(breaks = (-5:5)/100) +
	coord_cartesian(xlim = c(-0.05,1.05), ylim = c(-0.055,0.055), expand = F) +
	theme_few() +
	theme(aspect.ratio = 16/9,
				#legend.background = element_rect(fill = "grey90"),
				axis.text = element_text(size = 9),
				#panel.ontop = T, panel.background = element_rect(fill = NA),
				panel.grid = element_line(color = "grey80", size = 0.25, linetype = "solid"),
				panel.grid.major.x = element_line(color = "grey80", size = 0.25, linetype = "solid"),
				panel.grid.major.y = element_line(color = "grey80", size = 0.25, linetype = "solid"),
				panel.grid.minor = element_blank(),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text.x = element_text(face = "bold", size = 12),
				strip.text.y = element_text(face = "bold.italic", size = 14, angle = 0),
				# legend.position = "bottom",
				# legend.key.width = unit(2.5, "cm"),
				# legend.background = element_rect(fill = "grey90", colour = "black", size = 0.25),
				# legend.key.height = unit(0.5, "cm"),
				legend.position="none",
				legend.title = element_blank()) +
	ylab("Delta genotype pair frequency") +
	xlab("Allele frequency")
gg

ggsave(gg, filename = sprintf("_plot.hmm.emission.delta.%s.pdf", prefix), width = 8.5, height = 4.9)






cat("HMM emission probabilities ...\n")

tmp = as.data.frame(ibd)
tmp = split(tmp, cut.freq(AAF))
tmp = lapply(tmp, function(x) {
	x = colSums(x)
	as.data.table(as.list(x))
})
tag = names(tmp)
tmp = rbindlist(tmp)
rownames(tmp) = tag
ibd = tmp

tmp = as.data.frame(non)
tmp = split(tmp, cut.freq(AAF))
tmp = lapply(tmp, function(x) {
	x = colSums(x)
	as.data.table(as.list(x))
})
tag = names(tmp)
tmp = rbindlist(tmp)
rownames(tmp) = tag
non = tmp


frq = as.numeric(rownames(ibd))

emiss.ibd = NULL
for (p in names(ibd)) {
	tmp = approx(x = frq, y = ibd[[p]], xout = AAF, rule = 2)$y
	emiss.ibd = cbind(emiss.ibd, tmp)
}
colnames(emiss.ibd) = names(ibd)


frq = as.numeric(rownames(non))

emiss.non = NULL
for (p in names(non)) {
	tmp = approx(x = frq, y = non[[p]], xout = AAF, rule = 2)$y
	emiss.non = cbind(emiss.non, tmp)
}
colnames(emiss.non) = names(non)





identical(rownames(non), rownames(ibd))
identical(names(non), names(ibd))

rs.non = rowSums(non)
rs.ibd = rowSums(ibd)

emiss = data.table(Frequency = as.numeric(rownames(non)))

for (p in names(non)) {
	emiss = cbind(emiss, x = non[[p]] / rs.non)
}

for (p in names(ibd)) {
	emiss = cbind(emiss, x = ibd[[p]] / rs.ibd)
}

names(emiss) = c("Frequency",
								 "NON_00",
								 "NON_01",
								 "NON_02",
								 "NON_11",
								 "NON_12",
								 "NON_22",
								 "IBD_00",
								 "IBD_01",
								 "IBD_02",
								 "IBD_11",
								 "IBD_12",
								 "IBD_22")

for (col in names(emiss)){
	emiss[[col]] = sprintf("%.8f", emiss[[col]])
}

write.table(emiss, "HMM.emission.prob.txt", append = F, quote = F, sep = " ", row.names = F, col.names = T)



emiss.plot = NULL

rs = rowSums(ibd)
frq = as.numeric(rownames(ibd))
for (p in names(ibd)) {
	emiss.plot = rbind(emiss.plot, data.table(state = "IBD", freq = frq, pair = p, prob = ibd[[p]] / rs))
}

rs = rowSums(non)
frq = as.numeric(rownames(non))
for (p in names(non)) {
	emiss.plot = rbind(emiss.plot, data.table(state = "NON", freq = frq, pair = p, prob = non[[p]] / rs))
}



emp.emiss = array(0, c(2, 6, length(AAF)), dimnames = list(c("IBD", "NON"), names(ibd), NULL))

for (i in 1:length(AAF)) {
	emp.emiss["IBD", , i] = (emiss.ibd[i, ] + 1e-8) / sum(emiss.ibd[i, ])
	emp.emiss["NON", , i] = (emiss.non[i, ] + 1e-8) / sum(emiss.non[i, ])
}




cat("HMM initial probabilities ...\n")

# make new fk index
fklist = 2:(5000 - 2) # c(2:25, seq(30, 45, by = 5), seq(50, 90, by = 10), seq(100, 500, by=50))

inits = NULL

for (fk in fklist) {
	idx = which(AAC == fk)
	
	if (length(idx) < 100)
		next
	
	tmp = lapply(idx, function(x) {
		g = which(G[x, ] == 1)
		if (length(g) < 2) return(NULL) # others may be autozygous, hence higher count
		as.data.table(cbind(x, AAC[x], t(combn(g, 2))))
	})
	tmp = rbindlist(tmp)
	if (nrow(tmp) < 100) {
		next
	}
	names(tmp) = c("index", "fk", "g0", "g1")
	
	g0 = matrix(c(tmp$index, tmp$g0), ncol = 2, byrow = F)
	g1 = matrix(c(tmp$index, tmp$g1), ncol = 2, byrow = F)
	
	g0 = sim[ g0 ]
	g1 = sim[ g1 ]
	
	x = which(g0 == 1 & g1 == 1)
	
	ibd = length(x) / nrow(tmp)
	non = 1 - ibd
	
	out = data.table(fk = fk, frq = fk / 5000, ibd = ibd, non = non)
	
	print(out)
	
	inits = rbind(inits, out)
}


header = c("Frequency",
					 "NON",
					 "IBD")

cat(sprintf("%s\n", paste(header, collapse = " ")), file = "HMM.initial.prob.txt", append = F)

for (i in 1:nrow(inits))
{
	cat(sprintf("%.8f %.8f %.8f\n", inits$frq[i], inits$non[i], inits$ibd[i]), file = "HMM.initial.prob.txt", append = T)
}




cat("Saving ...\n")

save(emp.emiss, emp.init, file = sprintf("result.hmm_probs.%s.RData", prefix))




stop()


### emiss

library(ggplot2)
library(ggthemes)

p = (0:500) / 500
q = 1 - p

e = rbind(data.table(state = "IBD", pair = "00", freq = q, prob = p^3),
					data.table(state = "IBD", pair = "01", freq = q, prob = 2 * p^2 * q),
					data.table(state = "IBD", pair = "02", freq = q, prob = 0),
					data.table(state = "IBD", pair = "11", freq = q, prob = (p^2 * q) + (p * q^2)),
					data.table(state = "IBD", pair = "12", freq = q, prob = 2 * p * q^2),
					data.table(state = "IBD", pair = "22", freq = q, prob = q^3),
					data.table(state = " NON ", pair = "00", freq = q, prob = p^4),
					data.table(state = " NON ", pair = "01", freq = q, prob = 4 * p^3 * q),
					data.table(state = " NON ", pair = "02", freq = q, prob = 2 * p^2 * q^2),
					data.table(state = " NON ", pair = "11", freq = q, prob = 4 * p^2 * q^2),
					data.table(state = " NON ", pair = "12", freq = q, prob = 4 * p * q^3),
					data.table(state = " NON ", pair = "22", freq = q, prob = q^4))

emiss.plot$state[which(emiss.plot$state == "NON")] = " NON "

col = c("00" = "#6082E5",
				"01" = "#46B29D",
				"02" = "#B571E3",
				"11" = "#DBAC12",
				"12" = "#969130",
				"22" = "#DE473A")

gg = ggplot(emiss.plot) + 
	facet_grid(.~state) +
	#geom_bin2d(aes(x = freq, y = prob, fill = pair), bins = c(200,200), alpha=0.9) +
	geom_point(aes(x = freq, y = prob, colour = pair), size=1.25, shape = 16, alpha=0.75) +
	geom_line(data = e, aes(x = freq, y = prob, group = pair), colour = "white", alpha = 0.5, size = 1.5) +
	geom_line(data = e, aes(x = freq, y = prob, group = pair), colour = "black", alpha = 1.0, size = 0.5) +
	#geom_line(data = e, aes(x = freq, y = prob, color = pair), alpha = 0.95, size = 0.5) +
	scale_color_manual(values = col) +
	scale_fill_manual(values = col) +
	scale_x_continuous(expand = c(0.025, 0)) +
	scale_y_continuous(expand = c(0.025, 0)) +
	guides(colour = guide_legend(override.aes = list(size=5))) + 
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.margin.x = unit(x = 0.5, units = "cm"),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank(),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency") +
	ylab("Frequency")
gg
ggsave(gg, filename = sprintf("_plot.hmm.emission.%s.pdf", prefix), width = 9, height = 4.35)


# emiss expectations (general)

gg = ggplot(e) + 
	facet_grid(.~state) +
	geom_line(aes(x = freq, y = prob, color = pair), size = 1) +
	scale_color_manual(values = col) +
	scale_fill_manual(values = col) +
	scale_x_continuous(expand = c(0.025, 0), breaks = (0:4)/4, labels = as.character((0:4)/4)) +
	scale_y_continuous(expand = c(0.025, 0)) +
	guides(colour = guide_legend(override.aes = list(size=5))) + 
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.margin.x = unit(x = 0.5, units = "cm"),
				axis.title.x = element_text(margin=margin(10,0,0,0)),
				axis.title.y = element_text(margin=margin(0,10,0,0)),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank(),
				panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				strip.text = element_text(face = "bold")) +
	xlab("Allele frequency") +
	ylab("Genotype pair Frequency")
gg
ggsave(gg, filename = sprintf("../_plot.genpair-expected.pdf", prefix), width = 9, height = 4.35)






# initals

gg = ggplot(init.plot) + 
	geom_line(aes(x = fk, y = prob, colour = state), alpha = 0.5) +
	geom_point(aes(x = fk, y = prob, colour = state)) +
	scale_x_continuous(breaks = sort(unique(init.plot$fk)), labels = sort(unique(init.plot$fk))) +
	scale_colour_manual(values = c(IBD = "darkorange", NON = "royalblue2")) +
	theme_few() +
	theme(aspect.ratio = 1/4,
				legend.title = element_blank(),
				panel.grid = element_line(colour = "grey80"),
				panel.grid.major = element_line(colour = "grey80"),
				panel.grid.minor = element_blank()) +
	xlab("fk") +
	ylab("Probability")

ggsave(gg, filename = sprintf("_plot.hmm.initial.%s.pdf", prefix), width = 12, height = 3)






