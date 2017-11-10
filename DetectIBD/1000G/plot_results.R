


library(data.table)
library(ggplot2)
library(ggthemes)

load("data.1000GP_Phase3_chr20.RData")


load("result1.packs.1000GP_Phase3_chr20.RData")
res1 = res
load("result2.packs.1000GP_Phase3_chr20.RData")
res2 = res

k = intersect(names(res1), names(res2))
res = rbind(res1[, k, with=F], res2[, k, with=F])


panel = read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T, stringsAsFactors = F)
panel = as.data.table(panel)



ss = data.table(x0 = (2:25)-0.5, x1 = (2:25)+0.5, y0 = 1e-04, y1 = 50)
ss = ss[which((2:25) %% 2 != 0), ]

mb = c(seq(0.0001, 0.001, by = 0.0001),seq(0.001, 0.01, by = 0.001),seq(0.01, 0.1, by = 0.01),seq(0.1, 1, by = 0.1),seq(1, 10, by = 1),seq(10, 100, by = 10))


l = length(POS)


# physical
fgt = ((res$fgt.rhs.pos - res$fgt.lhs.pos) + 1) / 1e6
dgt = ((res$dgt.rhs.pos - res$dgt.lhs.pos) + 1) / 1e6
ghm = ((res$ghmm.rhs.pos - res$ghmm.lhs.pos) + 1) / 1e6
hhm = ((res$hhmm.rhs.pos - res$hhmm.lhs.pos) + 1) / 1e6
hhi = ((res$hhmmi.rhs.pos - res$hhmmi.lhs.pos) + 1) / 1e6


# genetic
fgt = ((DIST[res$fgt.rhs.idx] - DIST[res$fgt.lhs.idx]) + 1e-08)
dgt = ((DIST[res$dgt.rhs.idx] - DIST[res$dgt.lhs.idx]) + 1e-08)
ghm = ((DIST[res$ghmm.rhs.idx] - DIST[res$ghmm.lhs.idx]) + 1e-08)
hhm = ((DIST[res$hhmm.rhs.idx] - DIST[res$hhmm.lhs.idx]) + 1e-08)
hhi = ((DIST[res$hhmmi.rhs.idx] - DIST[res$hhmmi.lhs.idx]) + 1e-08)


a = which(res$fgt.lhs.idx != 1 & res$fgt.rhs.idx != l);    
b = which(res$dgt.lhs.idx != 1 & res$dgt.rhs.idx != l);     
c = which(res$ghmm.lhs.idx != 1 & res$ghmm.rhs.idx != l);  
d = which(res$hhmm.lhs.idx != 1 & res$hhmm.rhs.idx != l);   
e = which(res$hhmmi.lhs.idx != 1 & res$hhmmi.rhs.idx != l); 




d = rbind(data.table(idx = res$index[a], fk = res$fk[a], lhs = res$fgt.lhs.idx[a], rhs = res$fgt.rhs.idx[a], x = fgt[a], type = "(a) FGT"),
					data.table(idx = res$index[b], fk = res$fk[b], lhs = res$dgt.lhs.idx[b], rhs = res$dgt.rhs.idx[b], x = dgt[b], type = "(b) DGT"),
					data.table(idx = res$index[c], fk = res$fk[c], lhs = res$ghmm.lhs.idx[c], rhs = res$ghmm.rhs.idx[c], x = ghm[c], type = "(c) G-HMM"),
					data.table(idx = res$index[d], fk = res$fk[d], lhs = res$hhmm.lhs.idx[d], rhs = res$hhmm.rhs.idx[d], x = hhm[d], type = "(d) H-HMM"),
					data.table(idx = res$index[e], fk = res$fk[e], lhs = res$hhmmi.lhs.idx[e], rhs = res$hhmmi.rhs.idx[e], x = hhi[e], type = "(e) H-HMM*"))




q = split(d, list(d$type, d$fk))

q = lapply(q, function(x) {
	if (nrow(x) < 50) return(NULL)
	z = x
	w = fivenum(x$x)
	z$low = w[2]
	z$med = w[3]
	z$upp = w[4]
	z[1,]
})
q = rbindlist(q)



len = ggplot(q) + 
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=fk, ymin=low, ymax=upp, color=type), size = 1.5, position=position_dodge(width = 0.9)) +
	geom_point(aes(x=fk, y=med, group=type), size = 5, shape="-", position=position_dodge(width = 0.9)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), labels = format(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), scientific = F), minor_breaks = mb) +
	scale_x_continuous(breaks = 2:25) +
	coord_cartesian(ylim = c(0.0025, 3), xlim = c(1.5, 25.5), expand = F) +
	#coord_cartesian(ylim = c(1e-08, 2.1), xlim = c(1.5, 25.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = "top", #c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "horizontal",
				#legend.text = element_text(size = 9),
				legend.key.height = unit(0.1, "cm"),
				#legend.background = element_rect(fill = "white", colour = "grey50", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(colour = "grey80", size = 1/3),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Allele count, k") +
	#ylab("Physical length (Mb)") +
	ylab("Genetic length (cM)") +
	guides(color = guide_legend(override.aes = list(size = 5)))
len






pops = table(panel$super_pop)

q = NULL
for (pop in c("AFR", "AMR", "EUR", "SAS")) {
	tmp = rbind(data.table(idx = res$index[a], fk = res$fk[a], x = fgt[a], type = "(a) FGT"),
							data.table(idx = res$index[b], fk = res$fk[b], x = dgt[b], type = "(b) DGT"),
							data.table(idx = res$index[c], fk = res$fk[c], x = ghm[c], type = "(c) G-HMM"),
							data.table(idx = res$index[d], fk = res$fk[d], x = hhm[d], type = "(d) H-HMM"),
							data.table(idx = res$index[e], fk = res$fk[e], x = hhi[e], type = "(e) H-HMM*"))
	tmp$pop = pop
	tmp$frq = round(leg[[pop]][tmp$idx] * pops[pop] * 2)
	q = rbind(q, tmp)
}
d = q


del = which(d$frq < 2)
if (length(del) > 0) {
	d = d[-del,]
}



q = split(d, list(d$type, d$pop, d$frq))

q = lapply(q, function(x) {
	if (nrow(x) < 100) return(NULL)
	z = x
	w = fivenum(x$x)
	z$low = w[2]
	z$med = w[3]
	z$upp = w[4]
	z[1,]
})
q = rbindlist(q)



len = ggplot(q) + 
	facet_grid(pop~.) +
	geom_rect(data = ss, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1), fill = "grey80", alpha = 0.5) +
	geom_linerange(aes(x=frq, ymin=low, ymax=upp, color=type), size = 1.5, position=position_dodge(width = 0.9)) +
	geom_point(aes(x=frq, y=med, group=type), size = 5, shape="-", position=position_dodge(width = 0.9)) +
	scale_y_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), labels = format(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10), scientific = F), minor_breaks = mb) +
	scale_x_continuous(breaks = 2:25) +
	coord_cartesian(ylim = c(0.0025, 3), xlim = c(1.5, 25.5), expand = F) +
	#coord_cartesian(ylim = c(1e-08, 2.1), xlim = c(1.5, 25.5), expand = F) +
	theme_few() +
	theme(panel.border=element_rect(fill = NA, colour = "black", size = 0.5),
				legend.title=element_blank(),
				legend.position = "top", #c(0.99, 0.99), legend.justification = c(1,1),
				legend.direction = "horizontal",
				#legend.text = element_text(size = 9),
				legend.key.height = unit(0.1, "cm"),
				#legend.background = element_rect(fill = "white", colour = "grey50", size = 0.5),
				strip.text = element_text(face = "bold", hjust = 0),
				panel.grid = element_line(colour = "grey80", size = 1/3),
				panel.grid.major.x = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", size = 1/2),
				panel.grid.minor.y = element_line(colour = "grey80", size = 1/3)) +
	xlab("Population-specific allele count, k") +
	#ylab("Physical length (Mb)") +
	ylab("Genetic length (cM)") +
	guides(color = guide_legend(override.aes = list(size = 5)))
len


ggsave(len, filename = sprintf("_plot.length-phy.pdf"), width = 12, height = 9)
ggsave(len, filename = sprintf("_plot.length-gen.pdf"), width = 12, height = 9)

ggsave(len, filename = sprintf("_plot.length-phy.linear.pdf"), width = 12, height = 9)
ggsave(len, filename = sprintf("_plot.length-gen.linear.pdf"), width = 12, height = 9)







load("hhmmi-stats.packs.1000GP_Phase3_chr20.RData")

stat$fk = AAC[stat$index]

ggplot(stat) +
	facet_grid(fk~side) +
	geom_histogram(aes(iter-1), binwidth = 1, alpha = 0.5, color = "black") +
	geom_hline(yintercept = 0, size = 0.5, color = "black") +
	scale_x_continuous(breaks = 0:10) +
	coord_cartesian(xlim = c(-0.5, 10)) +
	theme_bw() +
	xlab("Number of corrections")

ggsave(filename = sprintf("_plot.hhmmi-stat.pdf"), width = 7, height = 9)










