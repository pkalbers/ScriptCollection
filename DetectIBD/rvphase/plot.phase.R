


library(data.table)
library(ggplot2)
library(ggthemes)

serr <- function(x) {
	sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
}


load("phase.local.ALL.RData")
load("phase.stack.all.RData")




###
### indeterminate and else per fk
###

sub = sprintf("%04d", unique(est.stack$this))

### segment

se = est.segment[sub]
st = tru.segment[sub]

se = lapply(se, as.data.table)
se = rbindlist(se)

st = lapply(st, as.data.table)
st = rbindlist(st)

se = split(se, se$fk)
st = split(st, st$fk)


### cluster

ce = est.cluster[sub]
ct = tru.cluster[sub]

ce = lapply(ce, as.data.table)
ce = rbindlist(ce)
ce = split(ce, ce$fk)

ct = lapply(ct, as.data.table)
ct = rbindlist(ct)
ct = split(ct, ct$fk)


### stack group

ge = est.stack
gt = tru.stack

ge = split(ge, ge$fk)

gt = split(gt, gt$fk)



get.hets = function(d) {
	rbind(data.table(fk = d$fk[1], type = "Indeterminate   " , mn = mean(d$hht / d$het, na.rm = T), se = serr(d$hht / d$het)),
				data.table(fk = d$fk[1], type = "Correct   " ,       mn = mean((d$inf-d$err) / d$het, na.rm = T), se = serr((d$inf-d$err) / d$het)),
				data.table(fk = d$fk[1], type = "Incorrect   " ,     mn = mean(d$err / d$het, na.rm = T), se = serr(d$err / d$het)))
}


pse = rbindlist(lapply(se, get.hets))
pce = rbindlist(lapply(ce, get.hets))
pge = rbindlist(lapply(ge, get.hets))

pst = rbindlist(lapply(st, get.hets))
pct = rbindlist(lapply(ct, get.hets))
pgt = rbindlist(lapply(gt, get.hets))


pse$mode = "Segment"
pce$mode = "Stacked"
pge$mode = "Combined"

pst$mode = "Segment"
pct$mode = "Stacked"
pgt$mode = "Combined"

pe = rbind(pse, pce, pge)
pt = rbind(pst, pct, pgt)

pe$mode = factor(pe$mode, levels = c("Segment", "Stacked", "Combined"), ordered = T)
pt$mode = factor(pt$mode, levels = c("Segment", "Stacked", "Combined"), ordered = T)


d = pt
d = pe

gg = ggplot(d) +
	facet_grid(.~mode) +
	geom_ribbon(aes(x=fk, ymin=mn-se, ymax=mn+se, fill = type), size = 0.5, alpha = 1/3) +
	geom_line(aes(x=fk, y=mn, colour = type)) +
	geom_point(aes(x=fk, y=mn, colour = type), size = 1.5, show.legend = T) +
	scale_x_continuous(expand = c(0.02, 0), breaks = c(2, 5, 10, 15, 20, 25)) +
	scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 1, by = 0.1), minor_breaks = seq(0, 1, by = 0.05), labels = seq(0, 1, by = 0.1)*100) +
	scale_fill_manual(values =   c("Indeterminate   " = "grey50", "Correct   " = "limegreen", "Incorrect   " = "brown2")) +
	scale_colour_manual(values = c("Indeterminate   " = "grey50", "Correct   " = "limegreen", "Incorrect   " = "brown2")) +
	coord_cartesian(ylim = c(0, 1)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		#legend.justification = c(1, 1), legend.position = c(1, 1)-0.01,
		#legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.key.size = unit(0.5, "cm"),
		legend.position = "bottom", 
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Mean proportion ± SE (%)")
gg

ggsave(gg, filename = "_plot.phase.fk_prop_truth.pdf", height = 4.5, width = 9)
ggsave(gg, filename = "_plot.phase.fk_prop_detected.pdf", height = 4.5, width = 9)




#
# error per fk
#

get.errs = function(d) {
	rbind(data.table(fk = d$fk[1], type = " Overall error   " ,         mn = mean(d$err / d$het, na.rm = T), se = serr(d$hht / d$het)),
				data.table(fk = d$fk[1], type = "Excluding singletons   " ,  mn = mean((d$err-d$err.sin) / d$het, na.rm = T), se = serr((d$err-d$err.sin) / d$het)),
				data.table(fk = d$fk[1], type = "Excluding singletons and other focal alleles   " , mn = mean((d$err-d$err.sin-d$err.foc) / d$het, na.rm = T), se = serr((d$err-d$err.sin-d$err.foc) / d$het)))
}

cols = c(" Overall error   " = "black", "Excluding singletons   " = "magenta2", "Excluding singletons and other focal alleles   " = "cyan3")

pse = rbindlist(lapply(se, get.errs))
pce = rbindlist(lapply(ce, get.errs))
pge = rbindlist(lapply(ge, get.errs))

pst = rbindlist(lapply(st, get.errs))
pct = rbindlist(lapply(ct, get.errs))
pgt = rbindlist(lapply(gt, get.errs))


pse$mode = "Segment"
pce$mode = "Consensus"
pge$mode = "Combined"

pst$mode = "Segment"
pct$mode = "Consensus"
pgt$mode = "Combined"

pe = rbind(pse, pce, pge)
pt = rbind(pst, pct, pgt)

pe$mode = factor(pe$mode, levels = c("Segment", "Consensus", "Combined"), ordered = T)
pt$mode = factor(pt$mode, levels = c("Segment", "Consensus", "Combined"), ordered = T)


d = pt
d = pe


pt$tag = "True IBD  "
pe$tag = "Detected IBD    "

d = rbind(pt, pe)
d$tag = factor(d$tag, levels = c("True IBD  ", "Detected IBD    "), ordered = T)


del = which(d$mode == "Combined")
if (length(del) > 0)
	d = d[-del, ]


gg = ggplot(d) +
	facet_grid(.~mode) +
	geom_ribbon(aes(x=fk, ymin=mn-se, ymax=mn+se, fill = type, group = interaction(type, tag)), size = 0.5, alpha = 1/3) +
	geom_line(aes(x=fk, y=mn, colour = type, group = interaction(type, tag))) +
	geom_point(aes(x=fk, y=mn, group = interaction(type, tag)), colour = "white", size = 1.25) +
	geom_point(aes(x=fk, y=mn, colour = type, shape = tag), size = 1.5) +
	#geom_point(aes(x=fk, y=mn, group = interaction(type, tag)), shape = 21, colour = "white", size = 1) +
	scale_x_continuous(expand = c(0.02, 0), breaks = c(2, 5, 10, 15, 20, 25)) +
	scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 1, by = 0.01), minor_breaks = seq(0, 1, by = 0.005), labels = seq(0, 1, by = 0.01)*100) +
	scale_fill_manual(values =   cols) +
	scale_colour_manual(values = cols, guide = guide_legend(override.aes = list(shape=32))) +
	scale_shape_manual(values = c(19, 21)) +
	coord_cartesian(ylim = c(0, 0.125)) +
	theme_few() +
	theme(aspect.ratio = 1,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		#legend.justification = c(1, 1), legend.position = c(1, 1)-0.01,
		#legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.key.size = unit(0.5, "cm"),
		legend.position = "bottom", 
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Mean proportion ± SE (%)")
gg

ggsave(gg, filename = "_plot.phase.fk_error.pdf", height = 5.5, width = 9)


#
# error site frq
#

get.frqs = function(d) {
	
	# filter out underestimated segments
	if (any(d$underest)) {
		del = which(d$underest)
		d = d[-del, ]
	}
	
	f = as.numeric(unlist(strsplit(d$err.frq, ",", T)))
	n = table(f)
	x = n / sum(n)
	data.table(fk = d$fk[1], frq = as.numeric(names(n))/5000, dn = x)
}

sub = c("2", "5", "15", "25")

pe = rbindlist(lapply(se[sub], get.frqs))
pt = rbindlist(lapply(st[sub], get.frqs))


d = pt
d = pe

d$fk = factor(d$fk)
d$dn.f = as.numeric(d$dn.f)

gg = ggplot(d) +
	#geom_line(aes(x=frq, y=dn.N, colour = fk), size = 1.5, show.legend = F) +
	geom_point(aes(x=dn.f, y=dn.N, colour = fk)) +
	#scale_x_log10(expand = c(0.02, 0), breaks = c((2:25), 50, 500, 5000)/5000 ) +
	scale_y_log10(expand = c(0.02, 0), 
								breaks = c(seq(0, 0.01, by = 0.0025), seq(0, 0.1, by = 0.025), seq(0, 1, by = 0.25)), 
								labels = c(seq(0, 0.01, by = 0.0025), seq(0, 0.1, by = 0.025), seq(0, 1, by = 0.25))*100) +
	coord_cartesian(xlim = c(0, 50), ylim = c(1e-4, 1)) +
	theme_few() +
	theme(#aspect.ratio = 1/2,
		panel.grid = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.y = element_line(colour = "grey85", size=1/3),
		panel.grid.major.y = element_line(colour = "grey85", size=1/3),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.title = element_blank(),
		#legend.justification = c(1, 1), legend.position = c(1, 1)-0.01,
		#legend.background = element_rect(fill = "white", colour = "black", size = 0.25),
		legend.key.size = unit(0.5, "cm"),
		legend.position = "bottom", 
		panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
		strip.text = element_text(face = "bold")) +
	xlab("Allele count, k") + ylab("Mean proportion ± SE (%)")
gg











#####
stop()
##### 


library(data.table)


files = dir(pattern = "^phase\\.overlap\\.([0-9]+)\\.est\\.RData$")

tmp = list()

for (file in files) {
	print(file)
	
	indv = sub("^phase\\.overlap\\.([0-9]+)\\.est\\.RData$", "\\1", file)
	
	load(file)
	
	if (length(est.overlap) > 100) {
		est.overlap = est.overlap[ sample(names(est.overlap), 100) ]
	}
	
	est.overlap = lapply(est.overlap, function(x) {
		x$dis.frq = NULL
		as.data.table(x)
	})
	est.overlap = rbindlist(est.overlap)
	
	tmp[[indv]] = est.overlap
}

est.overlap = tmp



files = dir(pattern = "^phase\\.overlap\\.([0-9]+)\\.tru\\.RData$")

tmp = list()

for (file in files) {
	print(file)
	
	indv = sub("^phase\\.overlap\\.([0-9]+)\\.tru\\.RData$", "\\1", file)
	
	load(file)
	
	if (length(tru.overlap) > 100) {
		tru.overlap = tru.overlap[ sample(names(tru.overlap), 100) ]
	}
	
	tru.overlap = lapply(tru.overlap, function(x) {
		x$dis.frq = NULL
		as.data.table(x)
	})
	tru.overlap = rbindlist(tru.overlap)
	
	tmp[[indv]] = tru.overlap
}

tru.overlap = tmp


save(est.overlap, tru.overlap, file = "phase.overlap.SUB.RData")

