

library(data.table)

load("~/Research/DetectIBD/data.OutOfAfricaHapMap20.H.RData")


files = dir(pattern = "^ibd.+ccf\\.txt$")

for (file in files) {
	cat(file, "\n")
	d = fread(file, header = T, stringsAsFactors = F)
	
	x = which(d$Clock == "C")
	d = d[x, ]
	
	d = d[, c("MarkerID", "Fk", "SampleID0", "Chr0", "SampleID1", "Chr1", "Shared", "SegmentLHS", "SegmentRHS")]
	
	d = unique(d)
	
	d$Act0 = d$Chr0
	d$Act1 = d$Chr1
	
	x = which(d$Shared == 1)
	
	c1  = matrix(c(d$MarkerID + 1, 1 + d$SampleID0 * 2 + d$Chr0), ncol = 2, byrow = F)
	c1a = matrix(c(d$MarkerID + 1, 1 + d$SampleID0 * 2 + abs((d$Chr0 + 1) - 2)), ncol = 2, byrow = F)
	c2  = matrix(c(d$MarkerID + 1, 1 + d$SampleID1 * 2 + d$Chr1), ncol = 2, byrow = F)
	c2a = matrix(c(d$MarkerID + 1, 1 + d$SampleID1 * 2 + abs((d$Chr1 + 1) - 2)), ncol = 2, byrow = F)
	
	h1  = H[c1]
	h1a = H[c1a]
	h2  = H[c2]
	h2a = H[c2a]
	
	
	q = NULL
	
	x = which(d$Shared == 1   &   h1 == 1   &   h2 == 1)
	if (length(x) > 0) q = rbind(q, d[x, ])
	
	x = which(d$Shared == 1   &   h1 == 0 & h1a == 1   &   h2 == 1)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act0 = abs((tmp$Chr0 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	x = which(d$Shared == 1   &   h1 == 1   &   h2 == 0 & h2a == 1)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act1 = abs((tmp$Chr1 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	x = which(d$Shared == 1   &   h1 == 0 & h1a == 1   &   h2 == 0 & h2a == 1)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act0 = abs((tmp$Chr0 + 1) - 2)
		tmp$Act1 = abs((tmp$Chr1 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	x = which(d$Shared == 0   &   h1 == 1   &   h2 == 0)
	if (length(x) > 0) q = rbind(q, d[x, ])
	
	x = which(d$Shared == 0   &   h1 == 0 & h1a == 1   &   h2 == 0)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act0 = abs((tmp$Chr0 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	x = which(d$Shared == 0   &   h1 == 1   &   h2 == 1 & h2a == 0)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act1 = abs((tmp$Chr1 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	x = which(d$Shared == 0   &   h1 == 0 & h1a == 1   &   h2 == 1 & h2a == 0)
	if (length(x) > 0) {
		tmp = d[x, ]
		tmp$Act1 = abs((tmp$Chr1 + 1) - 2)
		q = rbind(q, tmp)
	}
	
	cat(nrow(d), " >> ", nrow(q), "\n")
	
	q = q[order(q$MarkerID, abs((q$Shared + 1) - 2), q$SampleID0, q$SampleID1),]
	
	
	p = split(q, cut(1:nrow(q), breaks = 40))
	
	for (i in 1:length(p)) {
		z = p[[i]]
		write.table(z, file = sprintf("./packs/_pack%02d_%s", i, file), append = F, quote = F, row.names = F, col.names = T)
	}
}



###############


library(data.table)


files = dir(pattern = "^true.+")


d = NULL

for (file in files) {
	cat(file, "\n")
	
	tmp = fread(file, header = T, stringsAsFactors = F)
	
	tmp$data = sub("^true_pack([0-9]+)_ibd_con_dis\\.([a-z]+)_([A-Z])\\.n5000_id_([A-Z]+)\\.ccf\\.txt$", "\\2", file)
	tmp$phsd = sub("^true_pack([0-9]+)_ibd_con_dis\\.([a-z]+)_([A-Z])\\.n5000_id_([A-Z]+)\\.ccf\\.txt$", "\\3", file)
	tmp$near = sub("^true_pack([0-9]+)_ibd_con_dis\\.([a-z]+)_([A-Z])\\.n5000_id_([A-Z]+)\\.ccf\\.txt$", "\\4", file)
	
	d = rbind(d, tmp)
	
}

save(d, file = "../_true_ibd_con_dis.RData")



##############


library(data.table)
library(ggplot2)
library(ggthemes)

rmsle = function(a, b) {
	n = length(a)
	a = log10(a)
	b = log10(b)
	sqrt(sum((a - b)^2) / n)
}


load("_true_ibd_con_dis.RData")

table(d$data, d$phsd, d$near)


x = which(abs(d$LHS - d$TrueLHS) == 1)
if (length(x) > 0) d$TrueLHS[x] = d$LHS[x]

x = which(abs(d$RHS - d$TrueRHS) == 1)
if (length(x) > 0) d$TrueRHS[x] = d$RHS[x]


M = fread("../truH.marker.txt", header = T, stringsAsFactors = F)

d$pos = M$Position[d$MarkerID + 1]

d$lhs_pos = M$Position[d$LHS + 1]
d$rhs_pos = M$Position[d$RHS + 1]
d$true_lhs_pos = M$Position[d$TrueLHS + 1]
d$true_rhs_pos = M$Position[d$TrueRHS + 1]

d$lhs_gen = M$GenDist[d$LHS + 1]
d$rhs_gen = M$GenDist[d$RHS + 1]
d$true_lhs_gen = M$GenDist[d$TrueLHS + 1]
d$true_rhs_gen = M$GenDist[d$TrueRHS + 1]

d$lhs_dist = (d$pos - d$lhs_pos) + 1
d$rhs_dist = (d$rhs_pos - d$pos) + 1
d$true_lhs_dist = (d$pos - d$true_lhs_pos) + 1
d$true_rhs_dist = (d$true_rhs_pos - d$pos) + 1



# q = split(d, list(d$data, d$phsd, d$near))
# q = lapply(q, function(d) {
# 	d$key = sprintf("%d %d %d", d$MarkerID, d$SampleID0, d$SampleID1)
# 	d
# })
# key = Reduce(intersect, lapply(q, function(d) d$key))
# q = lapply(q, function(d, k) {
# 	x = which(d$key %in% k)
# 	d[x,]
# }, key)
# q = rbindlist(q)
# table(q$data, q$phsd, q$near)
# table(q$Shared)
# d = q


d$cord = "Concordant pairs"
x = which(d$Shared == 0)
d$cord[x] = "Discordant pairs"

d$type = "(a) Before error"
x = which(d$data == "err")
d$type[x] = "(b) After error"

d$mode = " Simulated haplotypes "
x = which(d$phsd == "P")
d$mode[x] = "Phased haplotypes"

d$method = "Nearest neighbours"
x = which(d$near == "RD")
d$method[x] = "Randomly selected"


dd = d
dd = dd[order(dd$Fk), ]



p = data.table(data = c(d$type, d$type), 
							 phsd = c(d$mode, d$mode), 
							 near = c(d$method, d$method), 
							 cord = c(d$cord, d$cord), 
							 est  = c(d$lhs_dist, d$rhs_dist), 
							 tru  = c(d$true_lhs_dist, d$true_rhs_dist))


q = p[which(p$data == "(a) Before error" & p$near == "Nearest neighbours"), ]
q$cord[which(q$cord == "Discordant pairs")] = "Discordant pairs (nearest neighbours)"
qq = p[which(p$data == "(a) Before error" & p$near == "Randomly selected"), ]
qq = qq[which(qq$cord == "Discordant pairs"), ]
qq$cord = "Discordant pairs (random selection)"
q = rbind(q, qq)

q = p[which(p$data == "(b) After error" & p$near == "Nearest neighbours"), ]
q$cord[which(q$cord == "Discordant pairs")] = "Discordant pairs (nearest neighbours)"
qq = p[which(p$data == "(b) After error" & p$near == "Randomly selected"), ]
qq = qq[which(qq$cord == "Discordant pairs"), ]
qq$cord = "Discordant pairs (random selection)"
q = rbind(q, qq)




a = by(q, list(q$phsd, q$cord), function(x) cor(x$est, x$tru, method = 'p')^2)
b = by(q, list(q$phsd, q$cord), function(x) cor(x$est, x$tru, method = 's'))
c = by(q, list(q$phsd, q$cord), function(x) rmsle(x$est, x$tru))
n = by(q, list(q$phsd, q$cord), nrow)

a = melt(array(a, dim = dim(a), dimnames = dimnames(a)));  names(a) = c("phsd", "cord", "val")
b = melt(array(b, dim = dim(b), dimnames = dimnames(b)));  names(b) = c("phsd", "cord", "val")
c = melt(array(c, dim = dim(c), dimnames = dimnames(c)));  names(c) = c("phsd", "cord", "val")
n = melt(array(n, dim = dim(n), dimnames = dimnames(n)));  names(n) = c("phsd", "cord", "val")

a$str = sprintf("paste('Pearson\\'s ', r^2 == '%s')", format(round(a$val, 3), trim = F, width = 5, digits = 3))
b$str = sprintf("paste('Spearman\\'s ', r[S] == '%s')", format(round(b$val, 3), trim = F, width = 5, digits = 3))
c$str = sprintf("'RMSLE' == '%s'", format(round(c$val, 3), trim = F, width = 5, digits = 3))

a$y = 55; a$x = 6 * 10^7
b$y = 15; b$x = 6 * 10^7
c$y = 5;  c$x = 6 * 10^7


tm = c(1:9, (1:9)*10, (1:9)*10^2, (1:9)*10^3, (1:9)*10^4, (1:9)*10^5, (1:9)*10^6, (1:9)*10^7, (1:9)*10^8)
tm = rbind(data.table(x0 = 0.1, x1 = 1.8, y0 = tm, y1 = tm), data.table(y0 = 0.1, y1 = 1.8, x0 = tm, x1 = tm))

lm = 10^(1:7)
lm = rbind(data.table(x0 = 0.1, x1 = 2.2, y0 = lm, y1 = lm), data.table(y0 = 0.1, y1 = 2.2, x0 = lm, x1 = lm))

gg = ggplot(q) +
	facet_grid(phsd~cord) +
	geom_abline(slope = 1, intercept = c(0,0), color = "black", alpha = 1/2, size = 1/2) +
	geom_raster(aes(tru, est, fill = ( ..density.. / max(..density..) ) * 100), stat = "bin2d", bins = 100) +
	geom_tile(aes(tru, est, fill = ( ..density.. / max(..density..) ) * 100), stat = "bin2d", bins = 100) +
	geom_text(data = a, aes(x, y, label = str), size = 3, parse = T, hjust = 1) +
	geom_text(data = b, aes(x, y, label = str), size = 3, parse = T, hjust = 1) +
	geom_text(data = c, aes(x, y, label = str), size = 3, parse = T, hjust = 1) +
	geom_segment(data = tm, aes(x = x0, xend = x1, y = y0, yend = y1), color = "grey50") +
	geom_segment(data = lm, aes(x = x0, xend = x1, y = y0, yend = y1), color = "grey50") +
	scale_x_log10(breaks = 10^(1:7), labels = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")) + # format(10^(1:7), trim = T, scientific = T, big.mark = ',')) +
	scale_y_log10(breaks = 10^(1:7), labels = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")) + # format(10^(1:7), trim = T, scientific = F, big.mark = ',')) +
	scale_fill_gradientn(colours = c("grey90", "gold1", "orange", "orangered", "purple2", "navy"), trans = "log10", na.value = "grey90", limits = c(0.1, 100), 
											 breaks = c(0.001, 0.1, 1,10,100), labels = c("0.001%", "< 0.1%", "1%","10%","100%")) +
	coord_cartesian(xlim = c(1.5, 10^7 * 9), ylim = c(1.5, 10^7 * 9), expand = F) +
	theme_few() +
	theme(aspect.ratio = 1,
				axis.ticks = element_line(colour = "grey50"),
				panel.border = element_rect(fill = NA, colour = "black", size = 1/2),
				strip.text.x = element_text(face = "bold", size = 12),
				strip.text.y = element_text(size = 13),
				plot.title = element_text(face = "bold", size = 13),
				legend.title = element_blank(),
				legend.text = element_text(size = 8),
				legend.key.width = unit(4, "mm"),
				legend.justification = c(0,1), legend.position = c(0.01,1),
				legend.background = element_blank()) +
	xlab("True breakpoint distance") + ylab("Inferred breakpoint distance") +
	ggtitle(q$data[1])
gg$theme$panel.spacing.x = unit(c(0.75, 0.25), "cm")
gg

ggsave(gg, filename = "___plot.ibd_dis_con.tru.pdf", width = 11, height = 7.5)
ggsave(gg, filename = "___plot.ibd_dis_con.err.pdf", width = 11, height = 7.5)

ggsave(gg, filename = "__plot.ibd_dis_con.A.pdf", width = 12, height = 8)
ggsave(gg, filename = "__plot.ibd_dis_con.B.pdf", width = 12, height = 8)






p = data.table(fk = d$Fk, 
							 data = d$type, 
							 phsd = d$mode, 
							 near = d$method, 
							 cord = d$cord, 
							 est = ((d$rhs_pos - d$lhs_pos) + 1) , 
							 tru = ((d$true_rhs_pos - d$true_lhs_pos) + 1) )

p = data.table(fk = d$Fk, 
							 data = d$type, 
							 phsd = d$mode, 
							 near = d$method, 
							 cord = d$cord, 
							 est = ((d$rhs_gen - d$lhs_gen) + 1e-8) , 
							 tru = ((d$true_rhs_gen - d$true_lhs_gen) + 1e-8) )

q = p[which(p$data == "(a) Before error" & p$near == "Nearest neighbours"), ]
q$cord[which(q$cord == "Discordant pairs")] = "Discordant pairs (nearest neighbours)"
qq = p[which(p$data == "(a) Before error" & p$near == "Randomly selected"), ]
qq = qq[which(qq$cord == "Discordant pairs"), ]
qq$cord = "Discordant pairs (random selection)"
q = rbind(q, qq)

q = p[which(p$data == "(b) After error" & p$near == "Nearest neighbours"), ]
q$cord[which(q$cord == "Discordant pairs")] = "Discordant pairs (nearest neighbours)"
qq = p[which(p$data == "(b) After error" & p$near == "Randomly selected"), ]
qq = qq[which(qq$cord == "Discordant pairs"), ]
qq$cord = "Discordant pairs (random selection)"
q = rbind(q, qq)



a = by(q, list(q$phsd, q$cord), function(x) cor(x$est, x$tru, method = 'p')^2)
b = by(q, list(q$phsd, q$cord), function(x) cor(x$est, x$tru, method = 's'))
c = by(q, list(q$phsd, q$cord), function(x) rmsle(x$est, x$tru))
n = by(q, list(q$phsd, q$cord), nrow)

a = melt(array(a, dim = dim(a), dimnames = dimnames(a)));  names(a) = c("phsd", "cord", "val")
b = melt(array(b, dim = dim(b), dimnames = dimnames(b)));  names(b) = c("phsd", "cord", "val")
c = melt(array(c, dim = dim(c), dimnames = dimnames(c)));  names(c) = c("phsd", "cord", "val")
n = melt(array(n, dim = dim(n), dimnames = dimnames(n)));  names(n) = c("phsd", "cord", "val")

a$str = sprintf("paste('Pearson\\'s ', r^2 == '%s')", format(round(a$val, 3), trim = F, width = 5, digits = 3))
b$str = sprintf("paste('Spearman\\'s ', r[S] == '%s')", format(round(b$val, 3), trim = F, width = 5, digits = 3))
c$str = sprintf("'RMSLE' == '%s'", format(round(c$val, 3), trim = F, width = 5, digits = 3))

a$y = 15; a$x = 48.5
b$y = 6.5; b$x = 48.5
c$y = 3.25;  c$x = 48.5


gr = NULL
for (k in 2:51) gr = rbind(gr, data.table(x=k-0.5, y = 10^(1:7)/1e6))

gg = ggplot(q) +
	facet_grid(phsd~cord) +
	geom_linerange(aes(fk, tru), stat = "summary", size = 2, alpha = 1/4, color = "blue", fun.ymin = function(z) {quantile(z,0.25)}, fun.ymax = function(z) {quantile(z,0.75)}) +
	geom_point(aes(fk,     tru), stat = "summary", size = 1.5, alpha = 1/4, color = "blue", fun.y = median, shape=19) +
	geom_point(aes(fk,     tru), stat = "summary", size = 1, alpha = 1/2, color = "white", fun.y = median, shape=20) +
	geom_vline(xintercept = (3:50)-0.5, alpha = 1/10, color = "black") +
	geom_point(data = gr, aes(x, y), alpha = 1/4, color = "black", shape = 95, size = 1.5) +
	geom_linerange(aes(fk, est), stat = "summary", size = 1/3, alpha = 1, color = "black", fun.ymin = function(z) {quantile(z,0.25)}, fun.ymax = function(z) {quantile(z,0.75)}) +
	geom_point(aes(fk,     est), stat = "summary", size = 1, alpha = 1, color = "black", fun.y = median, shape=19) +
	geom_point(aes(fk,     est), stat = "summary", size = 1, alpha = 1, color = "white", fun.y = median, shape=20) +
	geom_text(data = a, aes(x, y, label = str), size = 3.5, parse = T, hjust = 1) +
	geom_text(data = b, aes(x, y, label = str), size = 3.5, parse = T, hjust = 1) +
	geom_text(data = c, aes(x, y, label = str), size = 3.5, parse = T, hjust = 1) +
	coord_cartesian(xlim = c(1.5, 50.5), ylim = c(150, 10^7 * 3)/1e6, expand = F) +
	scale_x_continuous(breaks = c(2, 5, (1:5)*10)) +
	scale_y_log10(breaks = 10^(1:7)/1e6, labels = format(round(10^(1:7)/1e6, 3), scientific = F, drop0trailing = T)) + 
	#scale_y_log10(breaks = 10^(1:7), labels = c("10b", "100b", "1Kb", "10Kb", "100Kb", "1Mb", "10Mb")) + 
	theme_few() +
	theme(aspect.ratio = 1,
				axis.text = element_text(size = 9),
				panel.border=element_rect(fill = NA, colour = "black", size = 1/2),
				strip.text.x = element_text(face = "bold", size = 12),
				strip.text.y = element_text(size = 13),
				legend.justification=c(0,0), legend.position=c(1-0.985,1-0.985),
				legend.background = element_rect(fill = "white", colour = "grey80", size = 1/2),
				plot.title = element_text(face = "bold", size = 13),
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.margin = margin(0,1,1.5,1, "mm"),
				legend.text = element_text(size = 8),
				legend.title = element_blank()) +
	ylab("Genetic length (cM)") + 
	#ylab("Physical length") + 
	xlab("Focal allele count") +
	ggtitle(q$data[1])
gg$theme$panel.spacing.x = unit(c(0.75, 0.25), "cm")
gg

ggsave(gg, filename = "___plot.ibd_dis_con.tru_genlen.pdf", width = 11, height = 7.5)
ggsave(gg, filename = "___plot.ibd_dis_con.err_genlen.pdf", width = 11, height = 7.5)

ggsave(gg, filename = "___plot.ibd_dis_con.tru_phylen.pdf", width = 11, height = 7.5)
ggsave(gg, filename = "___plot.ibd_dis_con.err_phylen.pdf", width = 11, height = 7.5)

ggsave(gg, filename = "__plot.ibd_dis_con.A_len.pdf", width = 12, height = 5.75)
ggsave(gg, filename = "__plot.ibd_dis_con.B_len.pdf", width = 12, height = 5.75)







z = d[which(d$Shared==1 & d$data == "tru" & d$near == "NN"), ]
table(z$phsd)
z = split(z, z$phsd)

z = lapply(z, function(d) {
	d$key = sprintf("%d %d %d %d %d", d$MarkerID, d$SampleID0, d$Chr0, d$SampleID1, d$Chr1)
	d
})
key = Reduce(intersect, lapply(z, function(d) unique(d$key)))
z = lapply(z, function(d, k) {
	x = which(d$key %in% k)
	d[x,]
}, key)


mean(z$H$MarkerID - z$H$LHS)
mean(z$H$MarkerID - z$H$TrueLHS)

mean(z$P$MarkerID - z$P$LHS)
mean(z$P$MarkerID - z$P$TrueLHS)





