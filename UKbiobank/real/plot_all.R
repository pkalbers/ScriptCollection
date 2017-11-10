

library(data.table)
library(ggplot2)
library(ggthemes)


dhg.file = "real.---.chr7_107323983.ibd.distr.DHG.txt"
hmm.file = "real.---.chr7_107323983.age.pairs.HMM.Rec.HardBreaks.txt"

age.pairs = "real.---.chr7_107323983.age.pairs.HMM.Rec.HardBreaks.distr.txt"
age.sites = "real.---.chr7_107323983.age.sites.HMM.Rec.HardBreaks.distr.txt"


dhg.file = "real.1249839.chr7_44575563.ibd.distr.DHG.txt"
hmm.file = "real.1249839.chr7_44575563.age.pairs.HMM.Rec.txt"

age.pairs = "real.1249839.chr7_44575563.age.pairs.HMM.Rec.distr.txt"
age.sites = "real.1249839.chr7_44575563.age.sites.HMM.Rec.distr.txt"



out.txt = sub(".+chr([0-9]+)_([0-9]+)\\..+", "Chr \\1, Pos \\2", age.sites)
out.str = sub(".+chr([0-9]+)_([0-9]+)\\..+", "chr_\\1_pos_\\2", age.sites)

hard = "llksurf"
if (grepl("HardBreaks", age.sites)) hard = "hardbrk"



markers = fread(sprintf("ukb_chr%s_c.marker.txt", sub(".+chr([0-9]+)_([0-9]+)\\..+", "\\1", age.sites)), header = T, stringsAsFactors = F)

pos = markers$Position
names(pos) = as.character(markers$MarkerID)

gen = markers$GenDist
names(gen) = as.character(markers$MarkerID)



d.dhg = fread(dhg.file, header = T, stringsAsFactors = F)
d.dhg$pair = sprintf("%d + %d", d.dhg$SampleID0, d.dhg$SampleID1)
d.dhg$pos_FocalID = pos[as.character(d.dhg$FocalID)]
d.dhg$pos_BreakID = pos[as.character(d.dhg$BreakID)]



s = read.table(age.sites, header = T, stringsAsFactors = F)
p = read.table(age.pairs, header = T, stringsAsFactors = F)

s = s[, -1]
p = p[, -1]

times = as.numeric(sub("^t(.+)$", "\\1", names(s)))
names(times) = names(s)


d.hmm = fread(hmm.file, header = T, stringsAsFactors = F)


d.hmm$plength = pos[d.hmm$SegmentRHS] - pos[d.hmm$SegmentLHS]
d.hmm$glength = gen[d.hmm$SegmentRHS] - gen[d.hmm$SegmentLHS]

o = which(d.hmm$Shared == 0)

shr = d.hmm[-o,]
oth = d.hmm[o,]

1 / median(shr$glength)
1 / median(oth$glength)

1 / median(shr$plength)
1 / median(oth$plength)


###
q = fread("~/Research/rvage/africa/t1000.errH.age.pairs.HMM.Rec.txt", header = T, stringsAsFactors = F)
q = q[which(q$MarkerID == sample(q$MarkerID, 1)),]
w = fread("~/Research/rvage/africa/errH.marker.txt", header = T, stringsAsFactors = F)
pos = w$Position
names(pos) = as.character(w$MarkerID)
gen = w$GenDist
names(gen) = as.character(w$MarkerID)
q$plength = pos[q$SegmentRHS] - pos[q$SegmentLHS]
q$glength = gen[q$SegmentRHS] - gen[q$SegmentLHS]
o = which(q$Shared == 0)
shr = q[-o,]
oth = q[o,]

1 / median(shr$glength)
1 / median(oth$glength)

1 / median(shr$plength)
1 / median(oth$plength)

###


o = which(d.hmm$Shared == 0)
d.hmm = d.hmm[-o, ]
p.shr = p[-o, ]
p.out = p[o, ]




d.hmm$pair = sprintf("%d + %d", d.hmm$SampleID0, d.hmm$SampleID1)
d.hmm$pos_MarkerID = pos[as.character(d.hmm$MarkerID)]
d.hmm$pos_SegmentLHS = pos[as.character(d.hmm$SegmentLHS)]
d.hmm$pos_SegmentRHS = pos[as.character(d.hmm$SegmentRHS)]


if (length(unique(d.dhg$pair)) > 40) {
	z = intersect(unique(d.dhg$pair), unique(d.hmm$pair))
	z = sample(z, 40)
	d.dhg = d.dhg[which(d.dhg$pair %in% z), ]
	d.hmm = d.hmm[which(d.hmm$pair %in% z), ]
}







dhg = ggplot(d.dhg) + 
	facet_grid(pair~.) +
	geom_vline(aes(xintercept = pos_BreakID), colour = "grey40") +
	geom_point(data = data.frame(x=d.dhg$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=d.dhg$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(1e6, 300e6, by = 1e6), labels = sprintf("%d", seq(1e6, 300e6, by = 1e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = c(0.85, 1.15) * d.dhg$pos_FocalID[1]) +
	theme_few() + 
	theme(legend.position = "non",
				aspect.ratio = 1/75,
				legend.title = element_blank(),
				axis.line = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.text.y = element_blank(),
				strip.text.y = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") + ylab("") + ggtitle(out.txt)

ggsave(dhg, filename = sprintf("_plot.dhg.%s.png", out.str), width = 12, height = 12)




hmm = ggplot(d.hmm) + 
	facet_grid(pair~.) +
	geom_rect(aes(xmin = pos_SegmentLHS, xmax = pos_SegmentRHS, ymin = 0, ymax = 1), fill = "orange") +
	geom_point(data = data.frame(x=d.hmm$pos_MarkerID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=d.hmm$pos_MarkerID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(1e6, 300e6, by = 1e6), labels = sprintf("%d", seq(1e6, 300e6, by = 1e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = c(0.85, 1.15) * d.hmm$pos_MarkerID[1]) +
	theme_few() + 
	theme(aspect.ratio = 1/75,
				legend.position = "none",
				panel.background = element_rect(fill = "royalblue3"),
				legend.title = element_blank(),
				axis.line.x = element_blank(),
				axis.line.y = element_blank(),
				axis.text.y = element_blank(),
				strip.text.y = element_blank(),
				axis.ticks.y = element_blank()) +
	xlab("Physical position (Mb)") + ylab("") + ggtitle(out.txt)

ggsave(hmm, filename = sprintf("_plot.hmm.%s.png", out.str), width = 12, height = 12)




p.shr$pair = as.character(1:nrow(p.shr))
p.out$pair = as.character(1:nrow(p.out))

p.shr = melt(p.shr, id.vars = c("pair"))
p.out = melt(p.out, id.vars = c("pair"))

p.shr$variable = times[p.shr$variable] * 2 * 10000
p.out$variable = times[p.out$variable] * 2 * 10000

p.shr$type = " sharer "
p.out$type = "outgroup"

p = rbind(p.shr, p.out)

s = melt(s)
s$variable = times[s$variable] * 2 * 10000



gg = ggplot(p) +
	geom_line(aes(x = variable, y = value, colour = type, group = interaction(type, pair)), size = 0.5, alpha = 0.5) +
	geom_line(data = s, aes(x = variable, y = value), size = 1) +
	scale_color_manual(values = c("outgroup"="dodgerblue", " sharer "="orangered")) +
	scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = (c("1", "10", "100", "1000", "10000", "100000"))) +
	coord_cartesian(xlim = c(1, 100000)) +
	theme_few() +
	theme(aspect.ratio = 1/4,
				legend.position="none") +
	xlab("Time (generations") + ylab("CCF") + ggtitle(out.txt)

ggsave(gg, filename = sprintf("_plot.ccf_%s.%s.png", hard, out.str), width = 12, height = 6)



