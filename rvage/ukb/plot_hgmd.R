

library(data.table)
library(ggplot2)
library(ggthemes)


load("~/Research/UKbiobank/hgmd/hgmd_ukbb.RData")


####

POS = list()
GEN = list()
for (i in 1:22) {
	markers = fread(sprintf("../ukb_chr%d_c.marker.txt", i), header = T, stringsAsFactors = F)
	
	m.pos = markers$Position
	names(m.pos) = as.character(markers$MarkerID)
	
	m.gen = markers$GenDist
	names(m.gen) = as.character(markers$MarkerID)
	
	POS[[as.character(i)]] = m.pos
	GEN[[as.character(i)]] = m.gen
}


files = dir(pattern = ".+age\\.pairs\\.HMM\\.Rec\\.txt$")

d = list()

for (file in files) {
	
	ref = as.numeric(sub("^hgmd\\.([0-9]+)\\..+", "\\1", file))
	rsid = sub(".+_(rs[0-9]+)$", "\\1", hgmd$Label[which(hgmd$MarkerID == ref)[1]])
	
	if (!grepl("^rs", rsid)) next
	
	print(rsid)
	z = fread(file, header = T, stringsAsFactors = F)
	c=sub(".+\\.chr([0-9]+)_.+", "\\1", file)
	
	z$chr = as.numeric(c)
	z$rsid = rsid
	
	z$Chr0=NULL
	z$Chr1=NULL
	z$S_LHS=NULL
	z$S_RHS=NULL
	z$Missing=NULL
	
	mid = as.character(z$MarkerID)
	lhs = as.character(z$SegmentLHS)
	rhs = as.character(z$SegmentRHS)
	
	z$PhyPosFocal = POS[[c]][mid]
	z$GenPosFocal = GEN[[c]][mid]
	
	z$PhyPosLHS = POS[[c]][lhs]
	z$PhyPosRHS = POS[[c]][rhs]
	z$PhyLength = (z$PhyPosRHS - z$PhyPosLHS) + 1
	
	z$GenPosLHS = GEN[[c]][lhs]
	z$GenPosRHS = GEN[[c]][rhs]
	z$GenLength = (z$GenPosRHS - z$GenPosLHS) 

	
	lhs = as.character(z$SegmentLHS - 1)
	rhs = as.character(z$SegmentRHS + 1)
	
	z$BrkPhyPosLHS = POS[[c]][lhs]; if (any(is.na(z$BrkPhyPosLHS))) { na = which(is.na(z$BrkPhyPosLHS)); z$BrkPhyPosLHS[na] = POS[[c]][ as.character(z$SegmentLHS[na]) ] }
	z$BrkPhyPosRHS = POS[[c]][rhs]; if (any(is.na(z$BrkPhyPosRHS))) { na = which(is.na(z$BrkPhyPosRHS)); z$BrkPhyPosRHS[na] = POS[[c]][ as.character(z$SegmentRHS[na]) ] }
	z$BrkPhyLength = (z$PhyPosRHS - z$PhyPosLHS) + 1
	
	z$BrkGenPosLHS = GEN[[c]][lhs]; if (any(is.na(z$BrkGenPosLHS))) { na = which(is.na(z$BrkGenPosLHS)); z$BrkGenPosLHS[na] = GEN[[c]][ as.character(z$SegmentLHS[na]) ] }
	z$BrkGenPosRHS = GEN[[c]][rhs]; if (any(is.na(z$BrkGenPosRHS))) { na = which(is.na(z$BrkGenPosRHS)); z$BrkGenPosRHS[na] = GEN[[c]][ as.character(z$SegmentRHS[na]) ] }
	z$BrkGenLength = (z$BrkGenPosRHS - z$BrkGenPosLHS) 
	
	
	z$MarkerID=NULL
	z$SegmentLHS=NULL
	z$SegmentRHS=NULL
	
	z$Shared = as.logical(z$Shared)
	
	d[[rsid]] = z
}

d = lapply(d, function(x) if (length(table(x$Shared)) != 2) return(NULL) else return(x) )
dd=lapply(d, as.data.table)
dd=rbindlist(dd)
d = split(dd, dd$rsid)
d = lapply(d, function(x) if (min(table(x$Shared)) < 2) return(NULL) else return(x) )
dd=lapply(d, as.data.table)
dd=rbindlist(dd)
d = split(dd, dd$rsid)


z = lapply(d, function(x){
	a=which(x$Shared)
	b=which(!x$Shared)
	if (length(a) < 2) return(NULL)
	if (length(b) < 2) return(NULL)
	data.table(rsid=x$rsid[1],
						 phy.concord = median(x$BrkPhyLength[a] / 1e6),
						 phy.discord = median(x$BrkPhyLength[b] / 1e6), 
						 gen.concord = median(x$BrkGenLength[a]),
						 gen.discord = median(x$BrkGenLength[b]))
})
z = rbindlist(z)

gg=ggplot(z) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(x = phy.concord, y = phy.discord), alpha=0.9) +
	coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01, 100)) +
	scale_x_log10() +
	scale_y_log10() +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Concordant pairs") +
	ylab("Discordant pairs") +
	ggtitle("Median genetic length (Mb)")
gg
ggsave(gg, filename = "_plot.median_phylength.brk.png", width = 6, height = 6)

gg=ggplot(z) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(x = gen.concord, y = gen.discord), alpha=0.9) +
	coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01, 100)) +
	scale_x_log10() +
	scale_y_log10() +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Concordant pairs") +
	ylab("Discordant pairs") +
	ggtitle("Median physical length (Mb)")
gg
ggsave(gg, filename = "_plot.median_genlength.brk.png", width = 6, height = 6)


len = lapply(d, as.data.frame)
save(len, file="ukb_length.list.RData")
len = as.data.frame(rbindlist(d))
save(len, file="ukb_length.HMM.brk.RData")


t(as.data.table(lapply(d, function(x) sapply(split(x, x$Shared), function(y) as.data.table(median(y$GenLength))))))->q
plot(q[, 2], q[, 1], log="xy")
abline(0, 1)


use.hrdbrk = F


tags = as.numeric(sort(unique(sub("^hgmd\\.(.+)\\.chr.+", "\\1", dir(pattern = "^hgmd\\.(.+)\\.chr.+")))))

for (ref in tags) {
	
	chr = hgmd$Chromosome[which(hgmd$MarkerID == ref)[1]]
	pos = hgmd$Position[which(hgmd$MarkerID == ref)[1]]
	rsid = sub(".+_(rs[0-9]+)$", "\\1", hgmd$Label[which(hgmd$MarkerID == ref)[1]])
	fk = hgmd$AlleleCount1[which(hgmd$MarkerID == ref)[1]]
	
	if (fk > 15 | fk < 3) next
	
	if (rsid == "---") next
	
	
	prefix = sprintf("hgmd.%d.chr%d_%d", ref, chr, pos)
	
	
	hmm.file = sprintf("%s.age.pairs.HMM.Rec.txt", prefix)
	
	age.pairs = sprintf("%s.age.pairs.HMM.Rec.distr.txt", prefix)
	age.sites = sprintf("%s.age.sites.HMM.Rec.distr.txt", prefix)
	
	ccf.type = "llksurf"
	
	if (use.hrdbrk) {
		age.pairs = sprintf("%s.age.pairs.HMM.Rec.HardBreaks.distr.txt", prefix)
		age.sites = sprintf("%s.age.sites.HMM.Rec.HardBreaks.distr.txt", prefix)
		
		ccf.type = "hardbrk"
	}
	
	out.txt = sprintf("Chr %d, position %d (%s)", chr, pos, rsid)
	out.str = sprintf("chr_%d_pos_%d", chr, pos)
	
	
	if (!file.exists(hmm.file) || !file.exists(age.pairs) || !file.exists(age.sites)) next
	
	
	print(ref)
	
	m.pos = POS[[as.character(chr)]]
	
	
	
	s = as.data.frame(fread(age.sites, header = T, stringsAsFactors = F))
	p = as.data.frame(fread(age.pairs, header = T, stringsAsFactors = F))
	
	s = s[, -1]
	p = p[, -1]
	
	times = as.numeric(sub("^t(.+)$", "\\1", names(s)))
	names(times) = names(s)
	
	
	
	d.hmm = fread(hmm.file, header = T, stringsAsFactors = F)
	
	d.hmm$pair = sprintf("%d + %d", d.hmm$SampleID0, d.hmm$SampleID1)
	d.hmm$pos_MarkerID = m.pos[as.character(d.hmm$MarkerID)]
	d.hmm$pos_SegmentLHS = m.pos[as.character(d.hmm$SegmentLHS)]
	d.hmm$pos_SegmentRHS = m.pos[as.character(d.hmm$SegmentRHS)]
	
	o = which(d.hmm$Shared == 0)
	shr = d.hmm[-o, ]
	out = d.hmm[o, ]
	p.shr = p[-o, ]
	p.out = p[o, ]
	
	
	
	p.shr$pair = as.character(1:nrow(p.shr))
	p.out$pair = as.character(1:nrow(p.out))
	
	p.shr = melt(p.shr, id.vars = c("pair"))
	p.out = melt(p.out, id.vars = c("pair"))
	
	p.shr$variable = times[p.shr$variable] * 2 * 10000
	p.out$variable = times[p.out$variable] * 2 * 10000
	
	p.shr$type = "Concordant pair"
	p.out$type = "Discordant pair"
	
	p = rbind(p.shr, p.out)
	
	s = melt(s)
	s$variable = times[s$variable] * 2 * 10000
	
	
	
	ccf = ggplot(p) +
		geom_line(aes(x = variable, y = value, colour = type, group = interaction(type, pair)), size = 0.5, alpha = 0.5) +
		geom_line(data = s, aes(x = variable, y = value), size = 1) +
		scale_color_manual(values = c("Discordant pair"="dodgerblue", "Concordant pair"="orangered")) +
		scale_x_log10(expand = c(0.025, 0), breaks = c(1, 10, 100, 1000, 10000, 100000), labels = (c("1", "10", "100", "1000", "10000", "100000"))) +
		coord_cartesian(xlim = c(1, 100000)) +
		theme_few() +
		theme(aspect.ratio = 1/4,
					legend.title = element_blank(),
					legend.position="bottom") +
		xlab("Time (generations)") + ylab("CCF") + ggtitle(out.txt)
	
	ggsave(ccf, filename = sprintf("BRCA1_plot.%s.ccf_%s.png", out.str, ccf.type), width = 12, height = 6)
	
	
}





###


pair = list()

for (file in dir(pattern = "^hgmd\\..+\\.age\\.pairs\\.HMM\\.Rec\\.txt$")) {
	cat(file, " ")
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	if (nrow(tmp) > 0) {
		pos = POS[[sub(".+\\.chr([0-9]+).+", "\\1", file)]]
		gen = GEN[[sub(".+\\.chr([0-9]+).+", "\\1", file)]]
		
		tmp$pos_SegmentLHS = pos[as.character(tmp$SegmentLHS)]
		tmp$pos_SegmentRHS = pos[as.character(tmp$SegmentRHS)]
		tmp$phy_length = (tmp$pos_SegmentRHS - tmp$pos_SegmentLHS) + 1
		
		tmp$gen_SegmentLHS = gen[as.character(tmp$SegmentLHS)]
		tmp$gen_SegmentRHS = gen[as.character(tmp$SegmentRHS)]
		tmp$gen_length = (tmp$gen_SegmentRHS - tmp$gen_SegmentLHS)
		
		pair[[ sub("^hgmd\\.([0-9]+)\\.chr.+", "\\1", file) ]] = tmp
	}
}

z = lapply(pair, function(p) {
	if (nrow(p) < 2) return(NULL)
	
	s0 = which(p$Shared == 0)
	s1 = which(p$Shared == 1)
	
	if (length(s0) < 2 || length(s1) < 2) return(NULL)
	
	data.table(
		phy.concord = median(p$phy_length[s1] / 1e6),
		phy.discord = median(p$phy_length[s0] / 1e6),
		gen.concord = median(p$gen_length[s1]),
		gen.discord = median(p$gen_length[s0])
	)
})
z=rbindlist(z)


gg=ggplot(z) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(x = phy.concord, y = phy.discord)) +
	coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01, 100)) +
	scale_x_log10() +
	scale_y_log10() +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Concordant pairs") +
	ylab("Discordant pairs") +
	ggtitle("Median physical length (Mb)")
gg
ggsave(gg, filename = "_plot.median_phylength.png", width = 6, height = 6)

gg=ggplot(z) +
	geom_abline(slope = 1, intercept = 0) +
	geom_point(aes(x = gen.concord, y = gen.discord)) +
	coord_cartesian(xlim = c(0.005, 50), ylim = c(0.005, 50)) +
	scale_x_log10() +
	scale_y_log10() +
	theme_few() +
	theme(aspect.ratio = 1) +
	xlab("Concordant pairs") +
	ylab("Discordant pairs") +
	ggtitle("Median genetic length (cM)")
gg
ggsave(gg, filename = "_plot.median_genlength.png", width = 6, height = 6)


z = melt(z)
p = z[grep("^phy", z$variable), ]
p$value = p$value / 1e6
g = z[grep("^gen", z$variable), ]

p$variable[grep("concord")] = "Concordant"
p$variable[grep("discord")] = "Discordant"

g$variable[grep("concord")] = "Concordant"
g$variable[grep("discord")] = "Discordant"


ggplot(p) + geom_freqpoly(aes(value, colour=variable), position="dodge")+xlab("Physical length (Mb)")

ggplot(g) + geom_freqpoly(aes(value, colour=variable), position="dodge")+xlab("Genetic length (cM)")



#####


res = NULL

pair = list()

for (file in dir(pattern = "^hgmd\\..+\\.age\\.sites\\.HMM\\.Rec\\.txt$")) {
	print(file)
	tmp = read.table(file, header = T, stringsAsFactors = F)
	if (nrow(tmp) > 0) {
		
		res = rbind(res, tmp)
		
	}
}



######



x = match(res$MarkerID, hgmd$MarkerID)

llk = cbind(res, as.data.frame(hgmd)[x, ])
hrd = cbind(res, as.data.frame(hgmd)[x, ])

save(llk, hrd, file="ukbb_targets.RData")

ggplot(llk) + geom_point(aes(x=AlleleCount1, y=PostMean * 2 * 10000)) + scale_x_log10() + scale_y_log10()
ggplot(hrd) + geom_point(aes(x=AlleleCount1, y=PostMean * 2 * 10000)) + scale_x_log10() + scale_y_log10()

cor(llk$AlleleCount1, llk$PostMean)
cor(hrd$AlleleCount1, hrd$PostMean)

###



use.hrdbrk = F

ref = sample(unique(hgmd$MarkerID), 1)

chr = hgmd$Chromosome[which(hgmd$MarkerID == ref)[1]]
pos = hgmd$Position[which(hgmd$MarkerID == ref)[1]]
rsid = sub(".+_(rs[0-9]+)$", "\\1", hgmd$Label[which(hgmd$MarkerID == ref)[1]])
fk = hgmd$AlleleCount1[which(hgmd$MarkerID == ref)[1]]

print(fk)


prefix = sprintf("hgmd.%d.chr%d_%d", ref, chr, pos)


dhg.file = sprintf("%s.ibd.distr.DHG.txt", prefix)
hmm.file = sprintf("%s.age.pairs.HMM.Rec.txt", prefix)

age.pairs = sprintf("%s.age.pairs.HMM.Rec.distr.txt", prefix)
age.sites = sprintf("%s.age.sites.HMM.Rec.distr.txt", prefix)

ccf.type = "llksurf"

if (use.hrdbrk) {
	age.pairs = sprintf("%s.age.pairs.HMM.Rec.HardBreaks.distr.txt", prefix)
	age.sites = sprintf("%s.age.sites.HMM.Rec.HardBreaks.distr.txt", prefix)

	ccf.type = "hardbrk"
}

out.txt = sprintf("Chr %d, position %d (%s)", chr, pos, rsid)
out.str = sprintf("chr_%d_pos_%d", chr, pos)


if (!file.exists(dhg.file) || !file.exists(hmm.file) || !file.exists(age.pairs) || !file.exists(age.sites)) stop("Not present")




markers = fread(sprintf("../ukb_chr%d_c.marker.txt", chr), header = T, stringsAsFactors = F)

m.pos = markers$Position
names(m.pos) = as.character(markers$MarkerID)

m.gen = markers$GenDist
names(m.gen) = as.character(markers$MarkerID)









s = as.data.frame(fread(age.sites, header = T, stringsAsFactors = F))
p = as.data.frame(fread(age.pairs, header = T, stringsAsFactors = F))

s = s[, -1]
p = p[, -1]

times = as.numeric(sub("^t(.+)$", "\\1", names(s)))
names(times) = names(s)



d.hmm = fread(hmm.file, header = T, stringsAsFactors = F)

d.hmm$pair = sprintf("%d + %d", d.hmm$SampleID0, d.hmm$SampleID1)
d.hmm$pos_MarkerID = m.pos[as.character(d.hmm$MarkerID)]
d.hmm$pos_SegmentLHS = m.pos[as.character(d.hmm$SegmentLHS)]
d.hmm$pos_SegmentRHS = m.pos[as.character(d.hmm$SegmentRHS)]

o = which(d.hmm$Shared == 0)
shr = d.hmm[-o, ]
out = d.hmm[o, ]
p.shr = p[-o, ]
p.out = p[o, ]


shr$type = "Concordant"
out$type = "Discordant"

z = rbind(shr, out)
z = split(z, z$type)

tag = "Concordant"
tag = "Discordant"

for (tag in names(z)) {
	if (length(z[[tag]]$pair) > 15) z[[tag]] = z[[tag]][sample(1:nrow(z[[tag]]), 15),]
	
hmm = ggplot(z[[tag]]) + 
	facet_grid(pair~.) +
	geom_rect(aes(xmin = pos_SegmentLHS, xmax = pos_SegmentRHS, ymin = 0, ymax = 1), fill = "orange") +
	geom_point(data = data.frame(x=d.hmm$pos_MarkerID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
	geom_point(data = data.frame(x=d.hmm$pos_MarkerID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
	scale_x_continuous(breaks = seq(0, 300e6, by = 5e6), labels = sprintf("%d", seq(0, 300e6, by = 5e6) / 1e6)) +
	coord_cartesian(expand = F, xlim = range(m.pos)) + #xlim = c(0.8, 1.2) * d.hmm$pos_MarkerID[1]) +
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
hmm

ggsave(hmm, filename = sprintf("BRCA1_plot.%s.hmm.%s.png", out.str, tag), width = 12, height = 12)
}




p.shr$pair = as.character(1:nrow(p.shr))
p.out$pair = as.character(1:nrow(p.out))

p.shr = melt(p.shr, id.vars = c("pair"))
p.out = melt(p.out, id.vars = c("pair"))

p.shr$variable = times[p.shr$variable] * 2 * 10000
p.out$variable = times[p.out$variable] * 2 * 10000

p.shr$type = "Concordant pair"
p.out$type = "Discordant pair"

p = rbind(p.shr, p.out)

s = melt(s)
s$variable = times[s$variable] * 2 * 10000



ccf = ggplot(p) +
	geom_line(aes(x = variable, y = value, colour = type, group = interaction(type, pair)), size = 0.5, alpha = 0.5) +
	geom_line(data = s, aes(x = variable, y = value), size = 1) +
	scale_color_manual(values = c("Discordant pair"="dodgerblue", "Concordant pair"="orangered")) +
	scale_x_log10(expand = c(0.025, 0), breaks = c(1, 10, 100, 1000, 10000, 100000), labels = (c("1", "10", "100", "1000", "10000", "100000"))) +
	coord_cartesian(xlim = c(1, 100000)) +
	theme_few() +
	theme(aspect.ratio = 1/4,
				legend.title = element_blank(),
				legend.position="bottom") +
	xlab("Time (generations)") + ylab("CCF") + ggtitle(out.txt)
ccf

ggsave(ccf, filename = sprintf("_plot.%s.ccf_%s.png", out.str, ccf.type), width = 12, height = 6)








if (file.exists(dhg.file)) {
	
	d.dhg = fread(dhg.file, header = T, stringsAsFactors = F)
	d.dhg$pair = sprintf("%d + %d", d.dhg$SampleID0, d.dhg$SampleID1)
	d.dhg$pos_FocalID = m.pos[as.character(d.dhg$FocalID)]
	d.dhg$pos_BreakID = m.pos[as.character(d.dhg$BreakID)]
	
	
	dhg = ggplot(d.dhg) + 
		facet_grid(pair~.) +
		geom_vline(aes(xintercept = pos_BreakID), colour = "grey40") +
		geom_point(data = data.frame(x=d.dhg$pos_FocalID[1], y=0.5), aes(x, y), size = 2.5, colour = "black") +
		geom_point(data = data.frame(x=d.dhg$pos_FocalID[1], y=0.5), aes(x, y), size = 0.75, colour = "white") +
		scale_x_continuous(breaks = seq(0, 300e6, by = 5e6), labels = sprintf("%d", seq(0, 300e6, by = 5e6) / 1e6)) +
		coord_cartesian(expand = F, xlim = range(m.pos)) + #, xlim = c(0.8, 1.2) * d.dhg$pos_FocalID[1]) +
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
	dhg
	
	ggsave(dhg, filename = sprintf("_plot.%s.dhg.png", out.str), width = 12, height = 12)
	
}








S = read.table("~/Research/UKbiobank/ukb_chr17_c.sample.txt", skip = 1)


ref = hgmd$MarkerID[ grep("rs528724202", hgmd$Label) ]


chr = hgmd$Chromosome[which(hgmd$MarkerID == ref)[1]]
pos = hgmd$Position[which(hgmd$MarkerID == ref)[1]]
rsid = sub(".+_(rs[0-9]+)$", "\\1", hgmd$Label[which(hgmd$MarkerID == ref)[1]])
fk = hgmd$AlleleCount1[which(hgmd$MarkerID == ref)[1]]


prefix = sprintf("hgmd.%d.chr%d_%d", ref, chr, pos)

hmm.file = sprintf("%s.age.pairs.HMM.Rec.txt", prefix)

d = fread(hmm.file, header = T, stringsAsFactors = F)
d = d[which(d$Shared==1),]

s = sort(unique(c(d$SampleID0, d$SampleID1)))

ids = S$V2[match(s, S$V1)]
cat(sort(ids), sep="\n")
#cat(ids, sep="\n", file = sprintf("~/Research/UKbiobank/_%d.ids", ref))




z=NULL
for (file in dir(pattern = glob2rx("*.age.sites.HMM.Rec.txt"))) {
	
	print(file)
	tmp=read.table(file, header = T, stringsAsFactors = F)
	z=rbind(z, tmp)
	
}



















