#
# analyse error
#


###
make.error.matrix = function(tab) {
	m = matrix(0, nrow = 3, ncol = 3, dimnames = list(as.character(0:2), as.character(0:2)))
	
	is0 = which(tab$true.gt == 0)
	is1 = which(tab$true.gt == 1)
	is2 = which(tab$true.gt == 2)
	
	m["0", "0"] = length(which(tab$call.gt[is0] == 0))
	m["0", "1"] = length(which(tab$call.gt[is0] == 1))
	m["0", "2"] = length(which(tab$call.gt[is0] == 2))
	
	m["1", "0"] = length(which(tab$call.gt[is1] == 0))
	m["1", "1"] = length(which(tab$call.gt[is1] == 1))
	m["1", "2"] = length(which(tab$call.gt[is1] == 2))
	
	m["2", "0"] = length(which(tab$call.gt[is2] == 0))
	m["2", "1"] = length(which(tab$call.gt[is2] == 1))
	m["2", "2"] = length(which(tab$call.gt[is2] == 2))
	
	m
}


###
make.error.frqdep = function(tab, bins = 500, frq.size = 5008) {
	del = which(is.na(tab$frq.alt))
	if (length(del) > 0) {
		tab = tab[-del, ]
	}
	
	#brk = seq(min(tab$frq.alt), max(tab$frq.alt), length.out = bins)
	brk = seq(0, 1, length.out = bins+1)
	idx = cut(tab$frq.alt, brk, include.lowest = T)
	
	#brk = c(-1, unique(round((((0:250)/250)^exp(1)) * frq.size)))
	#idx = cut(round(tab$frq.alt * frq.size), brk, include.lowest = T)
	
	tab = split(tab, idx)
	#tab = split(tab, tab$frq.alt)
	
	subf = function(tab) {
		if (nrow(tab) == 0) {
			return(NULL)
		}
		is0 = which(tab$true.gt == 0)
		is1 = which(tab$true.gt == 1)
		is2 = which(tab$true.gt == 2)
		
		n0 = length(is0)
		n1 = length(is1)
		n2 = length(is2)
		
		alt = mean(tab$frq.alt)
		ref = 1 - alt
		
		data.table(ref = ref,
							 alt = alt,
							 num = nrow(tab),
							 
							 from.0.to.0 = if (n0 == 0) 0 else length(which(tab$call.gt[is0] == 0)),
							 from.0.to.1 = if (n0 == 0) 0 else length(which(tab$call.gt[is0] == 1)),
							 from.0.to.2 = if (n0 == 0) 0 else length(which(tab$call.gt[is0] == 2)),
							 
							 from.1.to.0 = if (n1 == 0) 0 else length(which(tab$call.gt[is1] == 0)),
							 from.1.to.1 = if (n1 == 0) 0 else length(which(tab$call.gt[is1] == 1)),
							 from.1.to.2 = if (n1 == 0) 0 else length(which(tab$call.gt[is1] == 2)),
							 
							 from.2.to.0 = if (n2 == 0) 0 else length(which(tab$call.gt[is2] == 0)),
							 from.2.to.1 = if (n2 == 0) 0 else length(which(tab$call.gt[is2] == 1)),
							 from.2.to.2 = if (n2 == 0) 0 else length(which(tab$call.gt[is2] == 2)))
	}
	
	tab = lapply(tab, subf)
	tab = rbindlist(tab)
	tab$frq = brk[-1] - diff(brk)/2
	tab
}


###
empirical.error = function(d, plotting = F) {
	pd = NULL
	
	n0 = d$from.0.to.0 + d$from.0.to.1 + d$from.0.to.2
	n1 = d$from.1.to.0 + d$from.1.to.1 + d$from.1.to.2
	n2 = d$from.2.to.0 + d$from.2.to.1 + d$from.2.to.2
	
	if (plotting) {
		
		pd = rbind(
			data.table(freq = d$frq, true.gt = "0", call.gt = "0", prop = d$from.0.to.0 / n0, true.n = n0),
			data.table(freq = d$frq, true.gt = "0", call.gt = "1", prop = d$from.0.to.1 / n0, true.n = n0),
			data.table(freq = d$frq, true.gt = "0", call.gt = "2", prop = d$from.0.to.2 / n0, true.n = n0),
			
			data.table(freq = d$frq, true.gt = "1", call.gt = "0", prop = d$from.1.to.0 / n1, true.n = n1),
			data.table(freq = d$frq, true.gt = "1", call.gt = "1", prop = d$from.1.to.1 / n1, true.n = n1),
			data.table(freq = d$frq, true.gt = "1", call.gt = "2", prop = d$from.1.to.2 / n1, true.n = n1),
			
			data.table(freq = d$frq, true.gt = "2", call.gt = "0", prop = d$from.2.to.0 / n2, true.n = n2),
			data.table(freq = d$frq, true.gt = "2", call.gt = "1", prop = d$from.2.to.1 / n2, true.n = n2),
			data.table(freq = d$frq, true.gt = "2", call.gt = "2", prop = d$from.2.to.2 / n2, true.n = n2)
		)
		del = which(is.na(pd$prop))
		if (length(del) > 0) {
			pd = pd[-del, ]
		}
		
	} else {
		pd = data.table(freq = d$alt, 
										p00 = d$from.0.to.0 / n0, p01 = d$from.0.to.1 / n0, p02 = d$from.0.to.2 / n0,
										p10 = d$from.1.to.0 / n1, p11 = d$from.1.to.1 / n1, p12 = d$from.1.to.2 / n1,
										p20 = d$from.2.to.0 / n2, p21 = d$from.2.to.1 / n2, p22 = d$from.2.to.2 / n2
		)
		del = c(which(is.na(pd$p00)), which(is.na(pd$p01)), which(is.na(pd$p02)),
						which(is.na(pd$p10)), which(is.na(pd$p11)), which(is.na(pd$p12)),
						which(is.na(pd$p20)), which(is.na(pd$p21)), which(is.na(pd$p22)))
		if (length(del) > 0) {
			del = unique(del)
			pd = pd[-del, ]
		}
	}
	
	pd
}


###
approx.error = function(d, af = (0:5000) / 5000, plotting = F) {
	pd = NULL
	
	a00 = approx(d$alt, d$from.0.to.0, xout = af, rule = 2)
	a01 = approx(d$alt, d$from.0.to.1, xout = af, rule = 2)
	a02 = approx(d$alt, d$from.0.to.2, xout = af, rule = 2)
	
	a10 = approx(d$alt, d$from.1.to.0, xout = af, rule = 2)
	a11 = approx(d$alt, d$from.1.to.1, xout = af, rule = 2)
	a12 = approx(d$alt, d$from.1.to.2, xout = af, rule = 2)
	
	a20 = approx(d$alt, d$from.2.to.0, xout = af, rule = 2)
	a21 = approx(d$alt, d$from.2.to.1, xout = af, rule = 2)
	a22 = approx(d$alt, d$from.2.to.2, xout = af, rule = 2)
	
	n0 = a00$y + a01$y + a02$y
	n1 = a10$y + a11$y + a12$y
	n2 = a20$y + a21$y + a22$y
	
	if (plotting) {
		
		pd = rbind(
			data.table(freq = af, true.gt = "0", call.gt = "0", prop = a00$y / n0),
			data.table(freq = af, true.gt = "0", call.gt = "1", prop = a01$y / n0),
			data.table(freq = af, true.gt = "0", call.gt = "2", prop = a02$y / n0),
			
			data.table(freq = af, true.gt = "1", call.gt = "0", prop = a10$y / n1),
			data.table(freq = af, true.gt = "1", call.gt = "1", prop = a11$y / n1),
			data.table(freq = af, true.gt = "1", call.gt = "2", prop = a12$y / n1),
			
			data.table(freq = af, true.gt = "2", call.gt = "0", prop = a20$y / n2),
			data.table(freq = af, true.gt = "2", call.gt = "1", prop = a21$y / n2),
			data.table(freq = af, true.gt = "2", call.gt = "2", prop = a22$y / n2)
		)
		del = which(is.na(pd$prop))
		if (length(del) > 0) {
			pd = pd[-del, ]
		}
		
	} else {
		
		pd = data.table(freq = af, 
										p00 = a00$y / n0, p01 = a01$y / n0, p02 = a02$y / n0,
										p10 = a10$y / n1, p11 = a11$y / n1, p12 = a12$y / n1,
										p20 = a20$y / n2, p21 = a21$y / n2, p22 = a22$y / n2
		)
		del = c(which(is.na(pd$p00)), which(is.na(pd$p01)), which(is.na(pd$p02)),
						which(is.na(pd$p10)), which(is.na(pd$p11)), which(is.na(pd$p12)),
						which(is.na(pd$p20)), which(is.na(pd$p21)), which(is.na(pd$p22)))
		if (length(del) > 0) {
			del = unique(del)
			pd = pd[-del, ]
		}
		
	}
	
	pd
}


#
#
#


library(ggplot2)
library(data.table)
library(ggthemes)


args = commandArgs(T)

pro.file = args[1]


cat("Loading data ... ")
load(pro.file)
cat("OK\n")
cat(sprintf(" # sites = %d\n", nrow(p)))

if (any(p$call.is.assumed)) {
	del = which(p$call.is.assumed)
	p = p[-del, ]
	cat(sprintf(" # sites after removing assumed called sites = %d\n", nrow(p)))
}
	
cat("Profiling errors ...")

error.matrix = make.error.matrix(p)

frq = make.error.frqdep(p)

error.table.empirical = empirical.error(frq)
error.table.approx    = approx.error(frq)

cat("OK\n")

#
#


cat("Saving ... ")
save(error.matrix, error.table.empirical, error.table.approx,
		 file = sprintf("result.%s", pro.file))
cat("OK\n")

# load(sprintf("result.%s", pro.file))

#
#


cat("Plotting ... ")

de = empirical.error(frq, plotting = T)
da = approx.error(frq, af = (0:5000)/5000, plotting = T)

em = error.matrix / rowSums(error.matrix) * 100
em = rbind(data.table(true.gt = "0", call.gt = "0", txt = sprintf("%.4f%%", em["0", "0"])),
					 data.table(true.gt = "0", call.gt = "1", txt = sprintf("%.4f%%", em["0", "1"])),
					 data.table(true.gt = "0", call.gt = "2", txt = sprintf("%.4f%%", em["0", "2"])),
					 
					 data.table(true.gt = "1", call.gt = "0", txt = sprintf("%.4f%%", em["1", "0"])),
					 data.table(true.gt = "1", call.gt = "1", txt = sprintf("%.4f%%", em["1", "1"])),
					 data.table(true.gt = "1", call.gt = "2", txt = sprintf("%.4f%%", em["1", "2"])),
					 
					 data.table(true.gt = "2", call.gt = "0", txt = sprintf("%.4f%%", em["2", "0"])),
					 data.table(true.gt = "2", call.gt = "1", txt = sprintf("%.4f%%", em["2", "1"])),
					 data.table(true.gt = "2", call.gt = "2", txt = sprintf("%.4f%%", em["2", "2"])))


de$true.gt = sprintf("True genotype is %s", de$true.gt)
de$call.gt = sprintf("Observed genotype is %s", de$call.gt)

da$true.gt = sprintf("True genotype is %s", da$true.gt)
da$call.gt = sprintf("Observed genotype is %s", da$call.gt)

em$true.gt = sprintf("True genotype is %s", em$true.gt)
em$call.gt = sprintf("Observed genotype is %s", em$call.gt)


gg = ggplot(data = da) + 
	facet_grid(call.gt~true.gt) + 
	geom_point(data = da, aes(x = freq, y = prop), size = 0.75, shape=1, colour="grey70") +
	geom_point(data = de, aes(x = freq, y = prop), size = 0.75, shape=1, colour="blue") + 
	geom_label(data = em, aes(x = 0.5, y = 0.5, label = txt), size = 3.5) +
	scale_x_continuous(expand = c(0.01, 0.01), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_y_continuous(expand = c(0.01, 0.01)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				panel.grid = element_line(colour = "grey80", linetype = "solid", size = 0.25),
				panel.grid.major = element_line(colour = "grey80", linetype = "solid", size = 0.25),
				panel.grid.minor = element_blank(),
				panel.border = element_rect(fill = NA, color = "black", size = 0.5),
				axis.text = element_text(size = 8),
				legend.position = "top") +
	xlab("Allele frequency") +
	ylab("Relative error proportion")

ggsave(filename = sprintf("_plot.error.%s.pdf", sub("^[^\\.]+\\.(.+)\\.RData$", "\\1", pro.file)), plot = gg, width = 10, height = 10)

cat("OK\n")





stop()
#### alternative plotting:


load("profile.1000g.weak.NA12878.RData")
p78 = p
load("profile.1000g.strict.NA12878.RData")
s78 = p
rm(p)

m78 = make.error.matrix(p78);   m78/rowSums(m78)*100
sm78 = make.error.matrix(s78);   sm78/rowSums(sm78)*100

f78 = make.error.frqdep(p78, bins = 200);  sum(f78$num)
sf78 = make.error.frqdep(s78, bins = 200);  sum(sf78$num)


d = empirical.error(f78, plotting = T)
sd = empirical.error(sf78, plotting = T)

a = as.data.frame(empirical.error(f78))
b = as.data.frame(f78)
for (i in 1:9) {
	cat(cor(b[, 3+i], b[, 3], method = 's'), "\n")
}


setwd("../affy6/")

load("profile.affy6.weak.NA12878.RData")
p78 = p
load("profile.affy6.weak.NA12877.RData")
p77 = p
rm(p)

load("profile.affy6.strict.NA12878.RData")
s78 = p
load("profile.affy6.strict.NA12877.RData")
s77 = p
rm(p)


setwd("../omni/")

load("profile.omni.weak.NA12878.RData")
p78 = p
load("profile.omni.weak.NA12877.RData")
p77 = p
rm(p)

load("profile.omni.strict.NA12878.RData")
s78 = p
load("profile.omni.strict.NA12877.RData")
s77 = p
rm(p)



m78 = make.error.matrix(p78);   m78/rowSums(m78)*100
m77 = make.error.matrix(p77);   m77/rowSums(m77)*100
mat = make.error.matrix(rbind(p78, p77));   mat/rowSums(mat)*100

sm78 = make.error.matrix(s78);   sm78/rowSums(sm78)*100
sm77 = make.error.matrix(s77);   sm77/rowSums(sm77)*100
smat = make.error.matrix(rbind(s78, s77));   smat/rowSums(smat)*100


f78 = make.error.frqdep(p78, bins = 200);  sum(f78$num)
f77 = make.error.frqdep(p77, bins = 200);  sum(f78$num)
frq = make.error.frqdep(rbind(p78, p77), bins = 200)

sf78 = make.error.frqdep(s78, bins = 200);  sum(sf78$num)
sf77 = make.error.frqdep(s77, bins = 200);  sum(sf78$num)
sfrq = make.error.frqdep(rbind(s78, s77), bins = 200)



d = empirical.error(frq, plotting = T)
sd = empirical.error(sfrq, plotting = T)





key = sprintf("%.5f %s %s", d$freq, d$true.gt, d$call.gt)
skey = sprintf("%.5f %s %s", sd$freq, sd$true.gt, sd$call.gt)
int = intersect(key, skey)
cor(d$prop[which(key %in% int)], sd$prop[which(skey %in% int)])^2

dd = cbind(d[which(key %in% int), ], sd[which(skey %in% int), ])
tmp = names(dd)
tmp[6:10] = sprintf("%s.s", tmp[6:10])
names(dd) = tmp

identical(dd$freq, dd$freq.s)
identical(dd$true.gt, dd$true.gt.s)
identical(dd$call.gt, dd$call.gt.s)



cs = by(d, list(d$true.gt, d$call.gt), function(d) cor(d$prop, d$true.n, method = "s"))
cs = array(cs, dim(cs), dimnames(cs))
cp = by(d, list(d$true.gt, d$call.gt), function(d) cor(d$prop, d$true.n, method = "p"))
cp = array(cp, dim(cp), dimnames(cp))

cs = melt(cs)
names(cs) = c("true.gt", "call.gt", "value")
cs$true.gt = sprintf("bold(t[%d])", cs$true.gt)
cs$call.gt = sprintf("bold(c[%d])", cs$call.gt)
cs$value = sprintf('rho[c] ~ "=" ~ %.2f', cs$value)

cp = melt(cp)
names(cp) = c("true.gt", "call.gt", "value")
cp$true.gt = sprintf("bold(t[%d])", cp$true.gt)
cp$call.gt = sprintf("bold(c[%d])", cp$call.gt)
cp$value = sprintf('r[c] ~ "=" ~ %.2f', cp$value)


scs = by(sd, list(sd$true.gt, sd$call.gt), function(d) cor(d$prop, d$true.n, method = "s"))
scs = array(scs, dim(scs), dimnames(scs))
scp = by(sd, list(sd$true.gt, sd$call.gt), function(d) cor(d$prop, d$true.n, method = "p"))
scp = array(scp, dim(scp), dimnames(scp))

scs = melt(scs)
names(scs) = c("true.gt", "call.gt", "value")
scs$true.gt = sprintf("bold(t[%d])", scs$true.gt)
scs$call.gt = sprintf("bold(c[%d])", scs$call.gt)
scs$value = sprintf('rho[s] ~ "=" ~ %.2f', scs$value)

scp = melt(scp)
names(scp) = c("true.gt", "call.gt", "value")
scp$true.gt = sprintf("bold(t[%d])", scp$true.gt)
scp$call.gt = sprintf("bold(c[%d])", scp$call.gt)
scp$value = sprintf('r[s] ~ "=" ~ %.2f', scp$value)



d$type = "weak"
sd$type = "strict"

dsd = rbind(d, sd)


pd = dsd
pd$true.gt = sprintf("bold(t[%s])", pd$true.gt)
pd$call.gt = sprintf("bold(c[%s])", pd$call.gt)

gg = ggplot(data = pd) + 
	facet_grid(call.gt~true.gt, labeller = label_parsed) + 
	#geom_point(aes(x = freq, y = prop, shape = type, colour = log10(true.n)), size = 0.75) +
	#geom_point(aes(x = freq, y = prop, shape = type), size = 0.75) +
	geom_point(aes(x = freq, y = prop, colour = true.gt, shape = type), size = 0.75) +
	#geom_linerange(aes(x = freq, ymin = prop, ymax = prop.s)) +
	#geom_point(aes(x = freq, y = -0.05, colour = log10(true.n))) +
	#geom_point(aes(x = freq, y = prop), size = 0.75) +
	geom_text(data = cs, aes(label = value, x = 0.5, y = 0.5), nudge_x = -0.125, nudge_y = +0.125, parse = T, size = 3) +
	geom_text(data = cp, aes(label = value, x = 0.5, y = 0.5), nudge_x = -0.125, nudge_y = -0.125, parse = T, size = 3) +
	geom_text(data = scs, aes(label = value, x = 0.5, y = 0.5), nudge_x = 0.125, nudge_y = +0.125, parse = T, size = 3) +
	geom_text(data = scp, aes(label = value, x = 0.5, y = 0.5), nudge_x = 0.125, nudge_y = -0.125, parse = T, size = 3) +
	scale_x_continuous(expand = c(0.075/3, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
	scale_y_continuous(expand = c(0.075, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "", "0.5", "", "1")) +
	scale_colour_manual(values = c("slateblue4", "springgreen4", "tomato4")) +
	#scale_color_gradient(low = "red", high = "green") + 
	scale_shape_manual(values = c(4, 1)) +
	theme_few() +
	theme(#aspect.ratio = 1/3,
				legend.title = element_blank(),
				legend.position = "none",
				strip.text.x = element_text(face = "bold", angle = 0),
				strip.text.y = element_text(face = "bold", angle = 0),
				panel.grid = element_line(colour = "grey85", linetype = "solid", size = 0.25),
				panel.grid.major = element_line(colour = "grey85", linetype = "solid", size = 0.25),
				panel.grid.minor = element_blank(),
				panel.border = element_rect(fill = NA, color = "black", size = 0.5),
				axis.text = element_text(size = 8),
				plot.margin = unit(c(1, 1, 1, 1)/2, "mm")) +
	xlab("Allele frequency") +
	ylab("Observed proportion")

gg

ggsave(gg, filename = "_plot.generrprop.1000g.pdf", width = 9, height = 4)
ggsave(gg, filename = "_plot.generrprop.affy6.pdf", width = 9, height = 4)
ggsave(gg, filename = "_plot.generrprop.omni.pdf", width = 9, height = 4)



ggplot(data = pd) + 
	facet_grid(.~true.gt) + 
	geom_line(aes(x = freq, y = true.n, colour = type)) +
	scale_y_log10() 






n = rbind(data.table(Individual = "NA12878", true.gt = 0, f = f78$frq, n = f78$from.0.to.0 + f78$from.0.to.1 + f78$from.0.to.2),
					data.table(Individual = "NA12878", true.gt = 1, f = f78$frq, n = f78$from.1.to.0 + f78$from.1.to.1 + f78$from.1.to.2),
					data.table(Individual = "NA12878", true.gt = 2, f = f78$frq, n = f78$from.2.to.0 + f78$from.2.to.1 + f78$from.2.to.2),
					data.table(Individual = "NA12877", true.gt = 0, f = f77$frq, n = f77$from.0.to.0 + f77$from.0.to.1 + f77$from.0.to.2),
					data.table(Individual = "NA12877", true.gt = 1, f = f77$frq, n = f77$from.1.to.0 + f77$from.1.to.1 + f77$from.1.to.2),
					data.table(Individual = "NA12877", true.gt = 2, f = f77$frq, n = f77$from.2.to.0 + f77$from.2.to.1 + f77$from.2.to.2))


gg = ggplot(n) +
	facet_grid(.~true.gt) + 
	#geom_bar(aes(x=factor(f), y=n, fill = Individual), stat = "identity") +
	geom_line(aes(x=(f), y=n, colour = Individual)) +
	scale_y_log10() +
	coord_cartesian(ylim = c(0.1, 10000))

gg



n = rbind(data.table(Individual = "NA12878", call.gt = 0, f = f78$frq, n = f78$from.0.to.0 + f78$from.1.to.0 + f78$from.2.to.0),
					data.table(Individual = "NA12878", call.gt = 1, f = f78$frq, n = f78$from.0.to.1 + f78$from.1.to.1 + f78$from.2.to.1),
					data.table(Individual = "NA12878", call.gt = 2, f = f78$frq, n = f78$from.0.to.2 + f78$from.1.to.2 + f78$from.2.to.2),
					data.table(Individual = "NA12877", call.gt = 0, f = f77$frq, n = f77$from.0.to.0 + f77$from.1.to.0 + f77$from.2.to.0),
					data.table(Individual = "NA12877", call.gt = 1, f = f77$frq, n = f77$from.0.to.1 + f77$from.1.to.1 + f77$from.2.to.1),
					data.table(Individual = "NA12877", call.gt = 2, f = f77$frq, n = f77$from.0.to.2 + f77$from.1.to.2 + f77$from.2.to.2))


gg = ggplot(n) +
	facet_grid(.~call.gt) + 
	#geom_bar(aes(x=factor(f), y=n, fill = Individual), stat = "identity") +
	geom_line(aes(x=(f), y=n, colour = Individual)) +
	scale_y_log10() +
	coord_cartesian(ylim = c(0.1, 10000))

gg




