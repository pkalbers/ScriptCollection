#
# analyse error
#


frq.size = 5008


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
make.error.frqdep = function(tab) {
	del = which(is.na(tab$frq.alt))
	if (length(del) > 0) {
		tab = tab[-del, ]
	}

	brk = c(-1, unique(round((((0:250)/250)^exp(1)) * frq.size)))
	idx = cut(round(tab$frq.alt * frq.size), brk, include.lowest = T)
	
	tab = split(tab, idx)
	
	subf = function(tab) {
		if (nrow(tab) == 0) {
			return(NULL)
		}
		
		is0 = which(tab$true.gt == 0)
		is1 = which(tab$true.gt == 1)
		is2 = which(tab$true.gt == 2)
		
		alt = (length(is2) + (length(is1) / 2)) / nrow(tab)
		ref = 1 - alt
		
		data.table(ref = ref,
							 alt = alt,
							 
							 from.0.to.0 = if (length(is0) == 0) 0 else length(which(tab$call.gt[is0] == 0)),
							 from.0.to.1 = if (length(is0) == 0) 0 else length(which(tab$call.gt[is0] == 1)),
							 from.0.to.2 = if (length(is0) == 0) 0 else length(which(tab$call.gt[is0] == 2)),
							 
							 from.1.to.0 = if (length(is1) == 0) 0 else length(which(tab$call.gt[is1] == 0)),
							 from.1.to.1 = if (length(is1) == 0) 0 else length(which(tab$call.gt[is1] == 1)),
							 from.1.to.2 = if (length(is1) == 0) 0 else length(which(tab$call.gt[is1] == 2)),
							 
							 from.2.to.0 = if (length(is2) == 0) 0 else length(which(tab$call.gt[is2] == 0)),
							 from.2.to.1 = if (length(is2) == 0) 0 else length(which(tab$call.gt[is2] == 1)),
							 from.2.to.2 = if (length(is2) == 0) 0 else length(which(tab$call.gt[is2] == 2)))
	}
	
	tab = lapply(tab, subf)
	rbindlist(tab)
}


###
subsample.profile = function(tab) {
	
	semi.random.subsample = function(i, hwp.count, assumed) {
		n = length(i)
		if (n <= 1) {
			return(i)
		}
		if (hwp.count == 0) { 
			return(c()) 
		}
		if (hwp.count == n) { 
			return(i) 
		}
		if (any(assumed)) {
			select = c(sample(which(!assumed)), sample(which(assumed)))
			select = select[1:hwp.count]
			return(i[select])
		}
		sample(i, hwp.count)
	}
	
	hwp.subsample = function(x) {
		if (is.null(x)) {
			return(NULL)
		}
		if (nrow(x) <= 1) {
			return(x)
		}
		
		i0 = which(x$true.gt == 0)
		i1 = which(x$true.gt == 1)
		i2 = which(x$true.gt == 2)
		
		N = nrow(x)
		
		ref = median(x$frq.ref)
		alt = 1 - ref
		
		obs.ratio = c(length(i0),
									length(i1),
									length(i2)) / N
		
		hwp.ratio = c(ref^2,
									2 * ref * alt,
									alt^2)
		
		hwp.count = as.integer( min(obs.ratio / hwp.ratio, na.rm = T) * N * hwp.ratio )
		
		i0 = semi.random.subsample(i0, hwp.count[1], x$call.is.assumed[i0])
		i1 = semi.random.subsample(i1, hwp.count[2], x$call.is.assumed[i1])
		i2 = semi.random.subsample(i2, hwp.count[3], x$call.is.assumed[i2])
		
		hwp = c(i0, i1, i2)
		#x$HWP[hwp] = TRUE
		
		if (length(hwp) == 0) {
			return(NULL)
		}
		
		as.data.table(x[hwp, ])
	}
	
	
	del = which(is.na(tab$frq.alt))
	if (length(del) > 0) {
		tab = tab[-del, ]
	}
	
	brk = c(-1, unique(round((((0:250)/250)^exp(1)) * frq.size)))
	idx = cut(round(tab$frq.alt * frq.size), brk, include.lowest = T)
	tab = split(tab, idx)
	
	n = length(tab)
	k = 0
	for (tag in names(tab)) {
		k = k + 1
		cat(sprintf(" %d of %d (%d)\n", k, n, nrow(tab[[tag]])))
		tab[[tag]] =  hwp.subsample(tab[[tag]])
	}
	
	tab = rbindlist(tab)
	
	tab
}


###
empirical.error = function(d, plotting = F) {
	
	n0 = d$from.0.to.0 + d$from.0.to.1 + d$from.0.to.2
	n1 = d$from.1.to.0 + d$from.1.to.1 + d$from.1.to.2
	n2 = d$from.2.to.0 + d$from.2.to.1 + d$from.2.to.2
	
	if (plotting) {
		
		pd = rbind(
			data.table(freq = d$alt, true.gt = "0", call.gt = "0", prop = d$from.0.to.0 / n0),
			data.table(freq = d$alt, true.gt = "0", call.gt = "1", prop = d$from.0.to.1 / n0),
			data.table(freq = d$alt, true.gt = "0", call.gt = "2", prop = d$from.0.to.2 / n0),
			
			data.table(freq = d$alt, true.gt = "1", call.gt = "0", prop = d$from.1.to.0 / n1),
			data.table(freq = d$alt, true.gt = "1", call.gt = "1", prop = d$from.1.to.1 / n1),
			data.table(freq = d$alt, true.gt = "1", call.gt = "2", prop = d$from.1.to.2 / n1),
			
			data.table(freq = d$alt, true.gt = "2", call.gt = "0", prop = d$from.2.to.0 / n2),
			data.table(freq = d$alt, true.gt = "2", call.gt = "1", prop = d$from.2.to.1 / n2),
			data.table(freq = d$alt, true.gt = "2", call.gt = "2", prop = d$from.2.to.2 / n2)
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


###
model.error = function(mat, gf0, gf1, gf2, plotting = F) {
	
	af = gf2 + (gf1 / 2)
	
	rs = rowSums(mat)
	
	e = mat / rs
	p = rs / sum(rs)
	
	e["0", ] = e["0", ] / p 
	e["1", ] = e["1", ] / p 
	e["2", ] = e["2", ] / p 
	
	e00 = e["0", "0"];  e01 = e["0", "1"];  e02 = e["0", "2"]
 	e10 = e["1", "0"];  e11 = e["1", "1"];  e12 = e["1", "2"]
 	e20 = e["2", "0"];  e21 = e["2", "1"];  e22 = e["2", "2"]
	
	s0 = (e00 * gf0) + (e01 * gf1) + (e02 * gf2)
	s1 = (e10 * gf0) + (e11 * gf1) + (e12 * gf2)
	s2 = (e20 * gf0) + (e21 * gf1) + (e22 * gf2)
	
	m00 = (e00 * gf0) / s0
	m01 = (e01 * gf1) / s0
	m02 = (e02 * gf2) / s0
	
	m10 = (e10 * gf0) / s1
	m11 = (e11 * gf1) / s1
	m12 = (e12 * gf2) / s1
	
	m20 = (e20 * gf0) / s2
	m21 = (e21 * gf1) / s2
	m22 = (e22 * gf2) / s2
	
	if (plotting) {
		pd = rbind(
			data.table(freq = af, true.gt = "0", call.gt = "0", prop = m00 ),
			data.table(freq = af, true.gt = "0", call.gt = "1", prop = m01 ),
			data.table(freq = af, true.gt = "0", call.gt = "2", prop = m02 ),
			
			data.table(freq = af, true.gt = "1", call.gt = "0", prop = m10 ),
			data.table(freq = af, true.gt = "1", call.gt = "1", prop = m11 ),
			data.table(freq = af, true.gt = "1", call.gt = "2", prop = m12 ),
			
			data.table(freq = af, true.gt = "2", call.gt = "0", prop = m20 ),
			data.table(freq = af, true.gt = "2", call.gt = "1", prop = m21 ),
			data.table(freq = af, true.gt = "2", call.gt = "2", prop = m22 )
		)
	} else {
		pd = data.table(freq = af, 
										p00 = m00, p01 = m01, p02 = m02,
										p10 = m10, p11 = m11, p12 = m12,
										p20 = m20, p21 = m21, p22 = m22
		)
	}
	
	pd
}


get.genotype.frq = function(d) {
	gc0 = d$from.0.to.0 + d$from.0.to.1 + d$from.0.to.2
	gc1 = d$from.1.to.0 + d$from.1.to.1 + d$from.1.to.2
	gc2 = d$from.2.to.0 + d$from.2.to.1 + d$from.2.to.2
	
	sum = gc0 + gc1 + gc2
	
	gf0 = gc0 / sum
	gf1 = gc1 / sum
	gf2 = gc2 / sum
	
	list("0" = gf0, "1" = gf1, "2" = gf2)
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


err.mat.all.all = make.error.matrix(p)
err.frq.all.all = make.error.frqdep(p)


cat("Subsampling profile\n")
q = subsample.profile(p)
cat("\nOK\n")
cat(sprintf(" # sites = %d\n", nrow(q)))


err.mat.all.sub = make.error.matrix(q)
err.frq.all.sub = make.error.frqdep(q)


cat("Removing assumed genotypes in experiment data ... ")
del = which(p$call.is.assumed)
if (length(del) > 0) {
	p = p [-del, ]
}
cat("OK\n")
cat(sprintf(" # sites = %d\n", nrow(p)))


err.mat.obs.all = make.error.matrix(p)
err.frq.obs.all = make.error.frqdep(p)


cat("Subsampling profile\n")
q = subsample.profile(p)
cat("\nOK\n")
cat(sprintf(" # sites = %d\n", nrow(q)))


err.mat.obs.sub = make.error.matrix(q)
err.frq.obs.sub = make.error.frqdep(q)


#
#


cat("Saving ... ")
save(err.mat.all.all, err.mat.all.sub, err.mat.obs.all, err.mat.obs.sub,
		 err.frq.all.all, err.frq.all.sub, err.frq.obs.all, err.frq.obs.sub,
		 empirical.error, approx.error, model.error,
		 file = sprintf("result.%s", pro.file))
cat("OK\n")


#
#


cat("Plotting ... ")

gf.all.all = get.genotype.frq(err.frq.all.all)
gf.all.sub = get.genotype.frq(err.frq.all.sub)
gf.obs.all = get.genotype.frq(err.frq.obs.all)
gf.obs.sub = get.genotype.frq(err.frq.obs.sub)


### compare models
d = rbind(cbind(model.error(err.mat.all.all, gf.all.all[["0"]], gf.all.all[["1"]], gf.all.all[["2"]], plotting = T), type = "All matched variants"),
					cbind(model.error(err.mat.all.sub, gf.all.sub[["0"]], gf.all.sub[["1"]], gf.all.sub[["2"]], plotting = T), type = "All matched variants, frequency subsampled"),
					cbind(model.error(err.mat.obs.all, gf.obs.all[["0"]], gf.obs.all[["1"]], gf.obs.all[["2"]], plotting = T), type = "Observed variants"),
					cbind(model.error(err.mat.obs.sub, gf.obs.sub[["0"]], gf.obs.sub[["1"]], gf.obs.sub[["2"]], plotting = T), type = "Observed variants, frequency subsampled"))

d$true.gt = sprintf("True genotype is %s", d$true.gt)
d$call.gt = sprintf("Observed genotype is %s", d$call.gt)

gg = ggplot(data = d) + 
	facet_grid(call.gt~true.gt) + 
	geom_line(aes(x = freq, y = prop, linetype = type, colour = type)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Allele frequency") +
	ylab("Relative error proportion")

ggsave(filename = sprintf("_plot.model.%s.pdf", sub("^[^\\.]+\\.(.+)\\.RData$", "\\1", pro.file)), plot = gg, width = 10, height = 10.5)


### compare approxs
d = rbind(cbind(approx.error(err.frq.all.all, af = (0:5000)/5000, plotting = T), type = "All matched variants"),
					cbind(approx.error(err.frq.all.sub, af = (0:5000)/5000, plotting = T), type = "All matched variants, frequency subsampled"),
					cbind(approx.error(err.frq.obs.all, af = (0:5000)/5000, plotting = T), type = "Observed variants"),
					cbind(approx.error(err.frq.obs.sub, af = (0:5000)/5000, plotting = T), type = "Observed variants, frequency subsampled"))

d$true.gt = sprintf("True genotype is %s", d$true.gt)
d$call.gt = sprintf("Observed genotype is %s", d$call.gt)

gg = ggplot(data = d) + 
	facet_grid(call.gt~true.gt) + 
	geom_line(aes(x = freq, y = prop, linetype = type, colour = type)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Allele frequency") +
	ylab("Relative error proportion")

ggsave(filename = sprintf("_plot.approx.%s.pdf", sub("^[^\\.]+\\.(.+)\\.RData$", "\\1", pro.file)), plot = gg, width = 10, height = 10.5)


### compare empiricals
d = rbind(cbind(empirical.error(err.frq.all.all, plotting = T), type = "All matched variants"),
					cbind(empirical.error(err.frq.all.sub, plotting = T), type = "All matched variants, frequency subsampled"),
					cbind(empirical.error(err.frq.obs.all, plotting = T), type = "Observed variants"),
					cbind(empirical.error(err.frq.obs.sub, plotting = T), type = "Observed variants, frequency subsampled"))

d$true.gt = sprintf("True genotype is %s", d$true.gt)
d$call.gt = sprintf("Observed genotype is %s", d$call.gt)

gg = ggplot(data = d) + 
	facet_grid(call.gt~true.gt) + 
	geom_line(aes(x = freq, y = prop, linetype = type, colour = type)) +
	theme_few() +
	theme(aspect.ratio = 1,
				legend.title = element_blank(),
				legend.position = "top") +
	xlab("Allele frequency") +
	ylab("Relative error proportion")

ggsave(filename = sprintf("_plot.empirical.%s.pdf", sub("^[^\\.]+\\.(.+)\\.RData$", "\\1", pro.file)), plot = gg, width = 10, height = 10.5)


###
l = list(all.all = list(mat = err.mat.all.all, frq = err.frq.all.all, gf = gf.all.all),
				 all.sub = list(mat = err.mat.all.sub, frq = err.frq.all.sub, gf = gf.all.sub),
				 obs.all = list(mat = err.mat.obs.all, frq = err.frq.obs.all, gf = gf.obs.all),
				 obs.sub = list(mat = err.mat.obs.sub, frq = err.frq.obs.sub, gf = gf.obs.sub))

for (tag in names(l)) {
	
	de = empirical.error(l[[tag]]$frq, plotting = T)
	da = approx.error(l[[tag]]$frq, af = (0:5000)/5000, plotting = T)
	dm = model.error(l[[tag]]$mat, l[[tag]]$gf[["0"]], l[[tag]]$gf[["1"]], l[[tag]]$gf[["2"]], plotting = T)
	
	em = l[[tag]]$mat / rowSums(l[[tag]]$mat) * 100
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
	
	dm$true.gt = sprintf("True genotype is %s", dm$true.gt)
	dm$call.gt = sprintf("Observed genotype is %s", dm$call.gt)
	
	em$true.gt = sprintf("True genotype is %s", em$true.gt)
	em$call.gt = sprintf("Observed genotype is %s", em$call.gt)
	
	gg = ggplot(data = da) + 
		facet_grid(call.gt~true.gt) + 
		geom_point(data = da, aes(x = freq, y = prop), size = 0.75, shape=16, colour="orangered") +
		geom_point(data = de, aes(x = freq, y = prop), size = 1.25, shape=16, colour="blue") + 
		#geom_line(data = dm, aes(x = freq, y = prop), colour = "black") +
		geom_label(data = em, aes(x = 0.5, y = 0.5, label = txt), size = 3.5) +
		theme_few() +
		theme(aspect.ratio = 1,
					legend.title = element_blank(),
					legend.position = "top") +
		xlab("Allele frequency") +
		ylab("Relative error proportion")
	
	ggsave(filename = sprintf("_plot.error.%s.%s.pdf", tag, sub("^[^\\.]+\\.(.+)\\.RData$", "\\1", pro.file)), plot = gg, width = 10, height = 10.5)
	
}

cat("OK\n")




