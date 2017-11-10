
library(data.table)
library(ggplot2)


# true allele frq
p = (0:5000)/5000
q = 1-p


# allele error rate
w = 0.01
r = 1 - w


# genotype error rate
e00 = r^2
e01 = 2 * w * r
e02 = w^2

e10 = w * r
e11 = ((w^2) + (r^2))
e12 = w * r

e20 = w^2
e21 = 2 * w * r
e22 = r^2


# genotype error rate
e = array(NA, c(4, 4), list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))

e["00", ] = c( (r^2) , w * r , w * r , (w^2) )
e["01", ] = c( w * r , (r^2) , (w^2) , w * r )
e["10", ] = c( w * r , (w^2) , (r^2) , w * r )
e["11", ] = c( (w^2) , w * r , w * r , (r^2) )

a00 = p^2
a01 = p * q
a10 = p * q
a11 = q^2

x00 = (a00 * e["00", "00"]) + (a01 * e["01", "00"]) + (a10 * e["10", "00"]) + (a11 * e["11", "00"])
x01 = (a00 * e["00", "01"]) + (a01 * e["01", "01"]) + (a10 * e["10", "01"]) + (a11 * e["11", "01"])
x10 = (a00 * e["00", "10"]) + (a01 * e["01", "10"]) + (a10 * e["10", "10"]) + (a11 * e["11", "10"])
x11 = (a00 * e["00", "11"]) + (a01 * e["01", "11"]) + (a10 * e["10", "11"]) + (a11 * e["11", "11"])

y00 = (x00 * e["00", "00"]) + (x00 * e["00", "01"]) + (x00 * e["00", "10"]) + (x00 * e["00", "11"])
y01 = (x01 * e["01", "00"]) + (x01 * e["01", "01"]) + (x01 * e["01", "10"]) + (x01 * e["01", "11"])
y10 = (x10 * e["10", "00"]) + (x10 * e["10", "01"]) + (x10 * e["10", "10"]) + (x10 * e["10", "11"])
y11 = (x11 * e["11", "00"]) + (x11 * e["11", "01"]) + (x11 * e["11", "10"]) + (x11 * e["11", "11"])


b00 = (sapply(e[, "00"], function(x) a00 * x))
b01 = (sapply(e[, "01"], function(x) a01 * x))
b10 = (sapply(e[, "10"], function(x) a10 * x))
b11 = (sapply(e[, "11"], function(x) a11 * x))


f = function(x, y) {
	get(sprintf("b%s", x))[, y]
}

non = function(x, y) {
	f = get(sprintf("b%s", x))
	if (y == "00") return( f[, "00"]^2 )
	if (y == "01") return( 4 * sqrt(f[, "00"])^3 * sqrt(f[, "11"]) )
	if (y == "02") return( 2 * f[, "00"] * f[, "11"] )
	if (y == "11") return( 4 * f[, "00"] * f[, "11"] )
	if (y == "12") return( 4 * sqrt(f[, "00"]) * sqrt(f[, "11"])^3 )
	if (y == "22") return( f[, "11"]^2 )
}

non = function(x, y) {
	if (y == "00") return( b00[, x]^2 )
	if (y == "01") return( (sqrt(b01[, x]) * sqrt(b01[, x])^3) + (sqrt(b10[, x]) * sqrt(b10[, x])^3) + (2 * sqrt(b11[, x]) * sqrt(b00[, x])^3) )
	if (y == "02") return( (2 * b00[, x]) * ( b11[, x]) )
	if (y == "11") return( (2 * b01[, x] * b01[, x]) + (2 * b10[, x] * b10[, x]) )
	if (y == "12") return( (sqrt(b01[, x]) * sqrt(b01[, x])^3) + (sqrt(b10[, x]) * sqrt(b10[, x])^3) + (2 * sqrt(b00[, x]) * sqrt(b11[, x])^3) )
	if (y == "22") return( b11[, x]^2 )
}

non = function(x, y) {
	f = get(sprintf("b%s", x))
	if (y == "00") return( f[, '00']^2 )
	if (y == "01") return( (sqrt(f[, '01']) * sqrt(f[, '01'])^3) + (sqrt(f[, '10']) * sqrt(f[, '10'])^3) + (2 * sqrt(f[, '11']) * sqrt(f[, '00'])^3) )
	if (y == "02") return( (2 * f[, '00']) * ( f[, '11']) )
	if (y == "11") return( (2 * f[, '01'] * f[, '01']) + (2 * f[, '10'] * f[, '10']) )
	if (y == "12") return( (sqrt(f[, '01']) * sqrt(f[, '01'])^3) + (sqrt(f[, '10']) * sqrt(f[, '10'])^3) + (2 * sqrt(f[, '00']) * sqrt(f[, '11'])^3) )
	if (y == "22") return( f[, '11']^2 )
}

p.non = function(xy) {
	x = sub("([0-2])([0-2])", "\\1", xy)
	y = sub("([0-2])([0-2])", "\\2", xy)
	
	r = function(xx) non(xx,'00') + non(xx,'01') + non(xx,'02') + non(xx,'11') + non(xx,'12') + non(xx,'22')
	
	z = rep(1, length(p))
	
	if (x == '0') z = z * r('00')
	if (x == '1') z = z * (r('01') + r('10'))
	if (x == '2') z = z * r('11')
	
	if (y == '0') z = z * r('00')
	if (y == '1') z = z * (r('01') + r('10'))
	if (y == '2') z = z * r('11')
	
	plot(z)
	
	lines(get(sprintf("n%s", xy)), col='red')
}

p.non = function(y) {
	
	plot( non(y, '00') + non(y, '01') + non(y, '02') + non(y, '11') + non(y, '12') + non(y, '22') )
	
	lines(get(sprintf("n%s", y)), col='red')
}
	

s00 = (b00[, "00"] + b01[, "00"] + b10[, "00"] + b11[, "00"])
s01 = (b00[, "01"] + b01[, "01"] + b10[, "01"] + b11[, "01"])
s10 = (b00[, "10"] + b01[, "10"] + b10[, "10"] + b11[, "10"])
s11 = (b00[, "11"] + b01[, "11"] + b10[, "11"] + b11[, "11"])

l00 = (b00[, "00"] + b00[, "01"] + b00[, "10"] + b00[, "11"])
l01 = (b01[, "00"] + b01[, "01"] + b01[, "10"] + b01[, "11"])
l10 = (b10[, "00"] + b10[, "01"] + b10[, "10"] + b10[, "11"])
l11 = (b11[, "00"] + b11[, "01"] + b11[, "10"] + b11[, "11"])


r00 = (              b01[, "00"] + b10[, "00"] + b11[, "00"])
r01 = (b00[, "01"]               + b10[, "01"] + b11[, "01"])
r10 = (b00[, "10"] + b01[, "10"]               + b11[, "10"])
r11 = (b00[, "11"] + b01[, "11"] + b10[, "11"]              )

l00 = (              b00[, "01"] + b00[, "10"] + b00[, "11"])
l01 = (b01[, "00"]               + b01[, "10"] + b01[, "11"])
l10 = (b10[, "00"] + b10[, "01"]               + b10[, "11"])
l11 = (b11[, "00"] + b11[, "01"] + b11[, "10"]              )


n00 = p^4
n01 = 4 * p^3 * q
n02 = 2 * p^2 * q^2
n11 = 4 * p^2 * q^2
n12 = 4 * p * q^3
n22 = q^4

i00 = p^3
i01 = 2 * p^2 * q
i02 = 0
i11 = p * q
i12 = 2 * p * q^2
i22 = q^3


plot(non('00','00') + non('00','01') + non('00','02') + non('00','11') + non('00','12') + non('00','22'))
lines(n00, col='red')

plot( (n00 + (n00 * s00^2)) + (n00 + (n00 * s00^2)))

plot( n01 + s00 - l00 )
plot( n02 + s00 - l00 )
plot( n11 + s00 - l00 )
plot( (1 - n12) + (r01 + r10 + r11) - ((n12) + (l01 + l10 + l11)) )
plot( n22 + s11 - l11 )


# true genotype frq
g0 = p^2
g1 = 2 * p * q
g2 = q^2


# observed allele frq (see below)
x = p * (1 - (2 * w)) + w
y = 1 - x


# observed genotype frq
x0 = (g0 * e00) + (g1 * e10) + (g2 * e20)  #  x^2 = (p^2 * (1-w)^2) + (2 * p * (1-p) * w * (1-w)) + ((1-p)^2 * w^2)
x1 = (g0 * e01) + (g1 * e11) + (g2 * e21)  #  2 * x * (1-x) = (p^2 * 2 * w * (1-w)) + (2 * p * (1-p) * ((w^2) + ((1-w)^2))) + ((1-p)^2 * 2 * w * (1-w))
x2 = (g0 * e02) + (g1 * e12) + (g2 * e22)  #  (1-x)^2 = (p^2 * w^2) + (2 * p * (1-p) * w * (1-w)) + ((1-p)^2 * (1-w)^2)


y0 = (g0 * f00) + (g1 * f10) + (g2 * f20) 
y1 = (g0 * f01) + (g1 * f11) + (g2 * f21) 
y2 = (g0 * f02) + (g1 * f12) + (g2 * f22)


p00 = (g0 * f00) / y0 
p01 = (g1 * f10) / y0
p02 = (g2 * f20) / y0

p10 = (g0 * f01) / y1
p11 = (g1 * f11) / y1
p12 = (g2 * f21) / y1

p20 = (g0 * f02) / y2
p21 = (g1 * f12) / y2
p22 = (g2 * f22) / y2



# probability of true --> observed genotype

p00 = (g0 * e00) / x0  #  (p^2 * (1-w)^2) / ((p^2 * (1-w)^2) + (2 * p * (1-p) * w * (1-w)) + ((1-p)^2 * w^2))  =  (p^2 * (1-w)^2) / (p * (1 - (2 * w)) + w)^2
p01 = (g1 * e10) / x0
p02 = (g2 * e20) / x0

p10 = (g0 * e01) / x1
p11 = (g1 * e11) / x1
p12 = (g2 * e21) / x1

p20 = (g0 * e02) / x2
p21 = (g1 * e12) / x2
p22 = (g2 * e22) / x2

range(p00 + p01 + p02)

d = rbind(data.table(f = q, p = p00, from = "from 0", to = "to 0"),
					data.table(f = q, p = p01, from = "from 0", to = "to 1"),
					data.table(f = q, p = p02, from = "from 0", to = "to 2"),
					data.table(f = q, p = p10, from = "from 1", to = "to 0"),
					data.table(f = q, p = p11, from = "from 1", to = "to 1"),
					data.table(f = q, p = p12, from = "from 1", to = "to 2"),
					data.table(f = q, p = p20, from = "from 2", to = "to 0"),
					data.table(f = q, p = p21, from = "from 2", to = "to 1"),
					data.table(f = q, p = p22, from = "from 2", to = "to 2"))

ggplot(d) + facet_grid(to~from) + geom_line(aes(x=f, y=p))




d = rbind(data.table(f = q, p = y0^2,        e = g0^2,         gt="00"),
					data.table(f = q, p = 2 * y0 * y1, e = 2 * g0 * g1,  gt="01"),
					data.table(f = q, p = 2 * y0 * y2, e = 2 * g0 * g2,  gt="02"),
					data.table(f = q, p = y1^2,        e = g1^2,         gt="11"),
					data.table(f = q, p = 2 * y1 * y2, e = 2 * g1 * g2,  gt="12"),
					data.table(f = q, p = y2^2,        e = g2^2,         gt="22"))


d = rbind(data.table(f = q, p = z00^2,        e = g0^2,         gt="00"),
					data.table(f = q, p = 2 * ((z00 * z01)+(z00 * z10)), e = 2 * g0 * g1,  gt="01"),
					data.table(f = q, p = 2 * z00 * z11, e = 2 * g0 * g2,  gt="02"),
					data.table(f = q, p = z01^2 + z10^2,        e = g1^2,         gt="11"),
					data.table(f = q, p = 2 * ((z11 * z01)+(z11 * z10)), e = 2 * g1 * g2,  gt="12"),
					data.table(f = q, p = z11^2,        e = g2^2,         gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")











n00 = p^4
n01 = 4 * p^3 * q
n02 = 2 * p^2 * q^2
n11 = 4 * p^2 * q^2
n12 = 4 * p * q^3
n22 = q^4

i00 = p^3
i01 = 2 * p^2 * q
i02 = 0
i11 = p * q
i12 = 2 * p * q^2
i22 = q^3


a00 = (n00 * p00^2)     + (n11 * p10^2)     + (n22 * p20^2)
a01 = (n00 * p00 * p01) + (n01 * p00 * p11) + (n02 * p00 * p21) + (n11 * p10 * p11) + (n12 * p10 * p21) + (n22 * p20 * p21) +
     	(n00 * p01 * p00) + (n01 * p11 * p00) + (n02 * p01 * p20) + (n11 * p11 * p10) + (n12 * p11 * p20) + (n22 * p21 * p20)

a01 = (n01 * p00 * p11) + (n01 * p01 * p10) + (2 * (  n00 * p00 * p01  +  n00 * p01 * p00 + 
																											n02 * p00 * p21  +  n02 * p01 * p20 + 
																											n11 * p10 * p11  +  n11 * p11 * p10 +
																											n12 * p10 * p21  +  n12 * p11 * p20 + 
																											n22 * p20 * p21  +  n22 * p21 * p20 ))

a02 = (n02 * p00 * p22) + (n02 * p02 * p20) + (2 * (  n00 * p00 * p02  +  
																									 		n01 * p00 * p12  +  n01 * p02 * p10 + 
																											n11 * p10 * p12  +  n11 * p12 * p10 +
																											n12 * p10 * p22  +  n12 * p12 * p20 + 
																											n22 * p20 * p22   ))

	(n11 * p10 * p11) + (n11 * p10 * p11) + (n22 * p20 * p21) + (n22 * p20 * p21) +
	(n01 * p00 * p11) + (n02 * p00 * p21) + (n12 * p10 * p21) + 
	(n01 * p01 * p10) + (n02 * p01 * p20) + (n12 * p11 * p20)


# true genotype pair frq: NON state
g00.non = g0 * g0      #  p^4
g01.non = 2 * g0 * g1  #  4 * p^3 * q
g02.non = 2 * g0 * g2  #  2 * p^2 * q^2
g11.non = g1 * g1      #  4 * p^2 * q^2
g12.non = 2 * g1 * g2  #  4 * p * q^3
g22.non = g2 * g2      #  q^4

# true genotype pair frq: IBD state
g00.ibd = g0 * p             #  p^3
g01.ibd = 2 * g0 * q         #  2 * p^2 * q
g02.ibd = 0                  #  0
g11.ibd = g0 * q + (p * g2)  #  (p^2 * q) + (p * q^2)
g12.ibd = 2 * p * g2         #  2 * p * q^2
g22.ibd = g2 * q             #  q^3


# observed genotype pair frq: NON state
x00.non = x0 * x0      #  p^4
x01.non = 2 * x0 * x1  #  4 * p^3 * q
x02.non = 2 * x0 * x2  #  2 * p^2 * q^2
x11.non = x1 * x1      #  4 * p^2 * q^2
x12.non = 2 * x1 * x2  #  4 * p * q^3
x22.non = x2 * x2      #  q^4

# observed genotype pair frq: IBD state
x00.ibd = x0 * x             #  p^3
x01.ibd = 2 * x0 * y         #  2 * p^2 * q
x02.ibd = 0                  #  0
x11.ibd = x0 * y + (x * x2)  #  (p^2 * q) + (p * q^2)
x12.ibd = 2 * x * x2         #  2 * p * q^2
x22.ibd = x2 * y             #  q^3


d = rbind(data.table(f = y, p = x0, e = g0, gt="00"),
					data.table(f = y, p = x1, e = g1, gt="01"),
					data.table(f = y, p = x2, e = g2, gt="02"),
					data.table(f = y, p = x0, e = g0, gt="11"),
					data.table(f = y, p = x1, e = g1, gt="12"),
					data.table(f = y, p = x2, e = g2, gt="22"))



z0 = (g0 * y00) + (g1 * y01) + (g2 * y02)
z1 = (g0 * y10) + (g1 * y11) + (g2 * y12)
z2 = (g0 * y20) + (g1 * y21) + (g2 * y22)

range(z0 + z1 + z2)

d = rbind(data.table(f = y, p = x0, e = g0, gt="0"),
					data.table(f = y, p = x1, e = g1, gt="1"),
					data.table(f = y, p = x2, e = g2, gt="2"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")







e = function(..., err = 0.01) {
	i = c(...)
	sum(err^i * (1 - err)^(4-i))
}

p = (0:5000)/5000
q = 1-p

g0 = p^2
g1 = 2 * p * q
g2 = q^2

n00 = p^4
n01 = 4 * p^3 * q
n02 = 2 * p^2 * q^2
n11 = 4 * p^2 * q^2
n12 = 4 * p * q^3
n22 = q^4

i00 = p^3
i01 = 2 * p^2 * q
i02 = 0
i11 = p * q
i12 = 2 * p * q^2
i22 = q^3

range(i00 + i01 + i02 + i11 + i12 + i22)


x = 0.01
P = (p - x) / (1 - 2*x)
P = ( -1 * sqrt(-1 * (2*x-1)^3 * sqrt(n00 * (2*x-1)^2)) + (4*x^3) - (4*x^2) + x ) / (2*x-1)^3
Q = 1 - P

((p^4) * (1-x)^4)+((4 * p^3 * (1-p)) * (x*(1-x)^3))+((2 * p^2 * (1-p)^2) * (x^2)*(1-x)^2)+((4 * p^2 * (1-p)^2) * (x^2)*(1-x)^2)+((4 * p * (1-p)^3) * (x^3)*(1-x))+(((1-p)^4) * x^4)

en00 = (n00 * e(0))  +  (n01 * e(1))  +  (n02 * e(2))  +  (n11 * e(2)) +  (n12 * e(3)) +  (n22 * e(4))
en01 = (n00 * e(1, 1))  +  (n01 * e(0, 2))  +  (n02 * e(1))  +  (n11 * e(1, 1)) +  (n12 * e(2, 2)) +  (n22 * e(3))
en02 = (n00 * e(2))  +  (n01 * e(1, 3))  +  (n02 * e(0))  +  (n11 * e(2)) +  (n12 * e(1, 3)) +  (n22 * e(2))
en11 = (n00 * e(2, 2, 2, 2))  +  (n01 * e(1, 1, 1, 1))  +  (n02 * e(2))  +  (n11 * e(0, 2, 2)) +  (n12 * e(1, 1, 1, 1)) +  (n22 * e(2, 2, 2, 2))
en12 = (n00 * e(3, 3))  +  (n01 * e(2, 2, 2))  +  (n02 * e(1))  +  (n11 * e(1, 3)) +  (n12 * e(0, 2)) +  (n22 * e(1, 1))
en22 = (n00 * e(4))  +  (n01 * e(3))  +  (n02 * e(2))  +  (n11 * e(2)) +  (n12 * e(1)) +  (n22 * e(0))

d = rbind(data.table(f = q, p = en00, e = n00, gt="00"),
					data.table(f = q, p = en01, e = n01, gt="01"),
					data.table(f = q, p = en02, e = n02, gt="02"),
					data.table(f = q, p = en11, e = n11, gt="11"),
					data.table(f = q, p = en12, e = n12, gt="12"),
					data.table(f = q, p = en22, e = n22, gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")



ei00 = (i00 * e(0))  +           (i01 * e(1, 1))  +        (i02 * e(2, 2))  +           (i11 * e(2)) +        (i12 * e(3)) +           (i22 * e(4))
ei01 = (i00 * e(1, 1))  +        (i01 * e(0, 2))  +        (i02 * e(1, 3))  +           (i11 * e(1, 1)) +     (i12 * e(2, 2)) +        (i22 * e(3, 3))
ei02 = (i00 * e(2))  +           (i01 * e(1, 3))  +        (i02 * e(0))  +              (i11 * e(2)) +        (i12 * e(1, 3)) +        (i22 * e(2))
ei11 = (i00 * e(2, 2, 2, 2))  +  (i01 * e(1, 1, 1, 1))  +  (i02 * e(2, 2, 2, 2))  +     (i11 * e(0, 2, 2)) +  (i12 * e(1, 1, 1, 1)) +  (i22 * e(2, 2, 2, 2))
ei12 = (i00 * e(3, 3))  +        (i01 * e(2, 2))  +        (i02 * e(1, 3))  +           (n11 * e(1, 3)) +     (i12 * e(0, 2)) +        (i22 * e(1, 1))
ei22 = (i00 * e(4))  +           (i01 * e(3, 3))  +        (i02 * e(2, 2))  +           (i11 * e(2)) +        (i12 * e(1)) +           (i22 * e(0))


ei00 = (i00 * e(0)*1)  +        (i01 * e(1, 1)*1)  +        (i02 * e(2, 2)*1)  +  (i11 * e(2)*1) +      (i12 * e(3, 3)*1) +        (i22 * e(4)*1)
ei01 = (i00 * e(1,1)*1)  +      (i01 * e(0, 2, 2, 2)*1)  +  (i02 * e(1, 3)*1)  +  (i11 * e(1)*1) +      (i12 * e(2, 2, 2)*1) +     (i22 * e(3, 3)*1)
ei02 = (i00 * e(2)*1)  +        (i01 * e(1, 3)*1)  +        (i02 * e(0, 4)*1)  +  (i11 * e(2)*1) +      (i12 * e(1, 3)*1) +        (i22 * e(2)*1)
ei11 = (i00 * e(2,2,2,2)*1)  +  (i01 * e(1, 1, 1, 1)*1)  +  (i02 * e(2, 2)*1)  +  (i11 * e(0,2,2)*1) +  (i12 * e(1, 1, 1, 1)*1) +  (i22 * e(2, 2, 2, 2)*1)
ei12 = (i00 * e(3,3)*1)  +      (i01 * e(2, 2, 2)*1)  +     (i02 * e(1, 3)*1)  +  (n11 * e(1, 3)*1) +   (i12 * e(0, 2, 2, 2)*1) +  (i22 * e(1, 1)*1)
ei22 = (i00 * e(4)*1)  +        (i01 * e(3, 3)*1)  +        (i02 * e(2, 2)*1)  +  (i11 * e(2)*1) +      (i12 * e(1, 1)*1) +        (i22 * e(0)*1)


ei00 = (i00 * e(0)*1)  +  (i01 * e(1, 1)*1)  +        (i02 * e(2, 2)*1)  +  (i11 * e(2)*1) +      (i12 * e(3, 3)*1) +        (i22 * e(4)*1)
ei01 = (i00 * e(1)*1)  +  (i01 * e(0, 2)*1)  +  (i01 * e(2, 2)*1)  +  (i02 * e(1, 3)*1)  +  (i11 * e(1)*1) +      (i12 * e(2, 2, 2)*1) +     (i22 * e(3)*1)
ei02 = (i00 * e(2)*1)  +  (i01 * e(1, 3)*1)  +        (i02 * e(0, 4)*1)  +  (i11 * e(2)*1) +      (i12 * e(1, 3)*1) +        (i22 * e(2)*1)
ei11 = (i00 * e(2)*1)  +  (i01 * e(1, 1, 1, 1)*1)  +  (i02 * e(2, 2)*1)  +  (i11 * e(0,2,2)*1) +  (i12 * e(1, 1, 1, 1)*1) +  (i22 * e(2)*1)
ei12 = (i00 * e(3)*1)  +  (i01 * e(2, 2, 2)*1)  +     (i02 * e(1, 3)*1)  +  (n11 * e(1, 3)*1) +   (i12 * e(0, 2, 2, 2)*1) +  (i22 * e(1)*1)
ei22 = (i00 * e(4)*1)  +  (i01 * e(3, 3)*1)  +        (i02 * e(2, 2)*1)  +  (i11 * e(2)*1) +      (i12 * e(1, 1)*1) +        (i22 * e(0)*1)

range(ei00 + ei01 + ei02 + ei11 + ei12 + ei22)

d = rbind(data.table(f = q, p = ei00, e = i00, gt="00"),
					data.table(f = q, p = ei01, e = i01, gt="01"),
					data.table(f = q, p = ei02, e = i02, gt="02"),
					data.table(f = q, p = ei11, e = i11, gt="11"),
					data.table(f = q, p = ei12, e = i12, gt="12"),
					data.table(f = q, p = ei22, e = i22, gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")






p = (0:5000)/5000
q = 1-p

w = 0.01
r = 1 - w

e00 = 99.9420
e01 = 0.0410
e02 = 0.0169

e10 = 0.5495
e11 = 99.2814 
e12 = 0.1690 

e20 = 0.0332 
e21 = 0.2277 
e22 = 99.7391 


e00 = e00 / (e00 + e01 + e02)
e01 = e01 / (e00 + e01 + e02)
e02 = e02 / (e00 + e01 + e02)

e10 = e10 / (e10 + e11 + e12)
e11 = e11 / (e10 + e11 + e12) 
e12 = e12 / (e10 + e11 + e12) 

e20 = e20 / (e20 + e21 + e22) 
e21 = e21 / (e20 + e21 + e22) 
e22 = e22 / (e20 + e21 + e22) 



w = 0.01
r = 1 - w

e00 = r^2
e01 = 2 * w * r
e02 = w^2

e10 = 1 * w * r
e11 = 1 * ((w^2) + (r^2))
e12 = 1 * w * r

e20 = w^2
e21 = 2 * w * r
e22 = r^2



g0 = p^2
g1 = 2 * p * q
g2 = q^2



y0 = (g0 * e00) + (g0 * e01) + (g0 * e02)
y1 = (g1 * e10) + (g1 * e11) + (g1 * e12)
y2 = (g2 * e20) + (g2 * e21) + (g2 * e22)

x0 = (g0 * e00) + (g1 * e10) + (g2 * e20)
x1 = (g0 * e01) + (g1 * e11) + (g2 * e21)
x2 = (g0 * e02) + (g1 * e12) + (g2 * e22)


p00 = (g0 * e00) / x0
p01 = (g0 * e01) / x0
p02 = (g0 * e02) / x0
p10 = (g1 * e10) / x1
p11 = (g1 * e11) / x1
p12 = (g1 * e12) / x1
p20 = (g2 * e20) / x2
p21 = (g2 * e21) / x2
p22 = (g2 * e22) / x2

plot( ((p^4 * (p00*p00)) + (4*p^3*q * (4*p00*p10)) + (2*p^2*q^2 * (2*p00*p20)) + (4*p^2*q^2 * (4*p10*p10)) + (4*p*q^3 * (4*p10*p20)) + (q^4 * (p20*p20))) )
lines(p^4, col="red")

plot( (g0 * p10) + (g1 * p11) + (g2 * p12) )
lines(g1, col="red")

plot( (g2 * p22) + (g1 * p21) + (g0 * p20) )
lines( (g2 * p22) - (g2 * p21) - (g2 * p20) , col="blue")
lines(g2, col="red")


y0 = g0 * ((exp(g0 * e00) / x0) + (exp(g1 * e10) / x1) + (exp(g2 * e20) / x2))
y1 = g1 * ((exp(g0 * e01) / x0) + (exp(g1 * e11) / x1) + (exp(g2 * e21) / x2))
y2 = g2 * ((exp(g0 * e02) / x0) + (exp(g1 * e12) / x1) + (exp(g2 * e22) / x2))



range(y0 + y1 + y2)

d = rbind(data.table(f = q, p = y0, e = g0, gt="0"),
					data.table(f = q, p = y1, e = g1, gt="1"),
					data.table(f = q, p = y2, e = g2, gt="2"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype proportion")



d = rbind(data.table(f = q, p = y0^2,        e = g0^2,         gt="00"),
					data.table(f = q, p = 2 * y0 * y1, e = 2 * g0 * g1,  gt="01"),
					data.table(f = q, p = 2 * y0 * y2, e = 2 * g0 * g2,  gt="02"),
					data.table(f = q, p = y1^2,        e = g1^2,         gt="11"),
					data.table(f = q, p = 2 * y1 * y2, e = 2 * g1 * g2,  gt="12"),
					data.table(f = q, p = y2^2,        e = g2^2,         gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")










w = 0.01
r = 1 - w

Err = array(NA, c(3, 3), dimnames = list(0:2, 0:2))

Err['0', ] = c( r^2        ,  2 * w * r            ,  w^2 )
Err['1', ] = c( 2 * w * r  ,  2 * ((w^2) + (r^2))  ,  2 * w * r )
Err['2', ] = c( w^2        ,  2 * w * r            ,  r^2 )

err = function(i, j) {
	i = as.character(i)
	j = as.character(j)
	Err[i, j]
}


p = (0:5000)/5000
q = 1-p

Gen = array(NA, c(3, length(p)), dimnames = list(NULL, 0:2))

Gen['0', ] = p^2
Gen['1', ] = 2 * p * q
Gen['2', ] = q^2

gen = function(i, j = NULL) {
	if (is.null(j)) {
		return(Gen[as.character(i), ])
	}
	Gen[as.character(i), j]
}






w = 0.01
r = 1 - w

ErrNonPair = array(NA, c(6, 6), dimnames = list(c("00", "01", "02", "11", "12", "22"), c("00", "01", "02", "11", "12", "22")))

ErrNonPair['00', ] = c( r^4            ,  4 * r^3 * w                    ,  2 * r^2 * w^2                  ,  4 * r^2 * w^2                      ,  4 * r * w^3                    ,  w^4            )
ErrNonPair['01', ] = c( 4 * r^3 * w    ,  4 * (r^4 + (3 * r^2 * w^2))    ,  (4 * r^3 * w) + (4 * r * w^3)  ,  (8 * r^3 * w) + (8 * r * w^3)      ,  4 * (w^4 + (3 * r^2 * w^2))    ,  4 * w^3 * r    )
ErrNonPair['02', ] = c( 2 * r^2 * w^2  ,  (4 * r^3 * w) + (4 * r * w^3)  ,  2 * (r^4 + w^4)                ,  8 * r^2 * w^2                      ,  (4 * r^3 * w) + (4 * r * w^3)  ,  2 * w^2 * r^2  )
ErrNonPair['11', ] = c( 4 * r^2 * w^2  ,  (8 * r^3 * w) + (8 * r * w^3)  ,  8 * r^2 * w^2                  ,  4 * (r^4 + w^4 + (2 * r^2 * w^2))  ,  (8 * r^3 * w) + (8 * r * w^3)  ,  4 * w^2 * r^2  )
ErrNonPair['12', ] = c( 4 * r * w^3    ,  4 * (w^4 + (3 * r^2 * w^2))    ,  (4 * r^3 * w) + (4 * r * w^3)  ,  (8 * r^3 * w) + (8 * r * w^3)      ,  4 * (r^4 + (3 * r^2 * w^2))    ,  4 * w * r^3    )
ErrNonPair['22', ] = c( w^4            ,  4 * w^3 * r                    ,  2 * w^2 * r^2                  ,  4 * w^2 * r^2                      ,  4 * w * r^3                    ,  r^4            )

err.non.pair = function(a0, a1, b0, b1) {
	pt = sprintf("%d%d", a0, a1)
	pf = sprintf("%d%d", b0, b1)
	ErrNonPair[pt, pf]
}


ErrIbdPair = array(NA, c(6, 6), dimnames = list(c("00", "01", "02", "11", "12", "22"), c("00", "01", "02", "11", "12", "22")))

ErrIbdPair['00', ] = c( r^3                    ,  2 * r^2 * w                    ,  0  ,  (r^2 * w) + (r * w^2)        ,  2 * r * w^2                    ,  w^3                    )
ErrIbdPair['01', ] = c( 2 * r^2 * w            ,  2 * (r^3 + (r * w^2))          ,  0  ,  2 * ((r^2 * w) + (r * w^2))  ,  2 * (w^3 + (r^2 * w))          ,  2 * r * w^2            )
ErrIbdPair['02', ] = c( 0                      ,  0                              ,  0  ,  0                            ,  0                              ,  0                      )
ErrIbdPair['11', ] = c( (r * w^2) + (r^2 * w)  ,  (2 * r^2 * w) + (2 * r * w^2)  ,  0  ,  (2 * r^3) + (2 * w^3)        ,  (2 * r * w^2) + (2 * r^2 * w)  ,  (r^2 * w) + (r * w^2)  )
ErrIbdPair['12', ] = c( 2 * r * w^2            ,  2 * (w^3 + (r^2 * w))          ,  0  ,  2 * ((r^2 * w) + (r * w^2))  ,  2 * (r^3 + (r * w^2))          ,  2 * r^2 * w            )
ErrIbdPair['22', ] = c( w^3                    ,  2 * w^2 * r                    ,  0  ,  (r * w^2) + (r^2 * w)        ,  2 * r^2 * w                    ,  r^3                    )

err.ibd.pair = function(a0, a1, b0, b1) {
	pt = sprintf("%d%d", a0, a1)
	pf = sprintf("%d%d", b0, b1)
	ErrIbdPair[pt, pf]
}


p = (0:5000)/5000
q = 1-p

NonPair = array(NA, c(6, length(p)), dimnames = list(c("00", "01", "02", "11", "12", "22"), NULL))

NonPair['00', ] = p^4
NonPair['01', ] = 4 * p^3 * q
NonPair['02', ] = 2 * p^2 * q^2
NonPair['11', ] = 4 * p^2 * q^2
NonPair['12', ] = 4 * p * q^3
NonPair['22', ] = q^4

gen.non.pair = function(p0, p1, i = NULL) {
	if (is.null(i)) {
		return(NonPair[sprintf("%d%d", p0, p1), ])
	}
	NonPair[sprintf("%d%d", p0, p1), i]
}


IbdPair = array(NA, c(6, length(p)), dimnames = list(c("00", "01", "02", "11", "12", "22"), NULL))

IbdPair['00', ] = p^3
IbdPair['01', ] = 2 * p^2 * q
IbdPair['02', ] = 0
IbdPair['11', ] = (p^2 * q) + (p * q^2)
IbdPair['12', ] = 2 * p * q^2
IbdPair['22', ] = q^3

gen.ibd.pair = function(p0, p1, i = NULL) {
	if (is.null(i)) {
		return(IbdPair[sprintf("%d%d", p0, p1), ])
	}
	IbdPair[sprintf("%d%d", p0, p1), i]
}



g00 = gen.non.pair(0, 0)
g01 = gen.non.pair(0, 1)
g02 = gen.non.pair(0, 2)
g11 = gen.non.pair(1, 1)
g12 = gen.non.pair(1, 2)
g22 = gen.non.pair(2, 2)

s00 = 
	(g00 * err.non.pair(0, 0, 0, 0)) + 
	(g01 * err.non.pair(0, 1, 0, 0)) + 
	(g02 * err.non.pair(0, 2, 0, 0)) + 
	(g11 * err.non.pair(1, 1, 0, 0)) + 
	(g12 * err.non.pair(1, 2, 0, 0)) + 
	(g22 * err.non.pair(2, 2, 0, 0))

s01 = 
	(g00 * err.non.pair(0, 0, 0, 1)) + 
	(g01 * err.non.pair(0, 1, 0, 1)) + 
	(g02 * err.non.pair(0, 2, 0, 1)) + 
	(g11 * err.non.pair(1, 1, 0, 1)) + 
	(g12 * err.non.pair(1, 2, 0, 1)) + 
	(g22 * err.non.pair(2, 2, 0, 1))

s02 = 
	(g00 * err.non.pair(0, 0, 0, 2)) + 
	(g01 * err.non.pair(0, 1, 0, 2)) + 
	(g02 * err.non.pair(0, 2, 0, 2)) + 
	(g11 * err.non.pair(1, 1, 0, 2)) + 
	(g12 * err.non.pair(1, 2, 0, 2)) + 
	(g22 * err.non.pair(2, 2, 0, 2))

s11 = 
	(g00 * err.non.pair(0, 0, 1, 1)) + 
	(g01 * err.non.pair(0, 1, 1, 1)) + 
	(g02 * err.non.pair(0, 2, 1, 1)) + 
	(g11 * err.non.pair(1, 1, 1, 1)) + 
	(g12 * err.non.pair(1, 2, 1, 1)) + 
	(g22 * err.non.pair(2, 2, 1, 1))

s12 = 
	(g00 * err.non.pair(0, 0, 1, 2)) + 
	(g01 * err.non.pair(0, 1, 1, 2)) + 
	(g02 * err.non.pair(0, 2, 1, 2)) + 
	(g11 * err.non.pair(1, 1, 1, 2)) + 
	(g12 * err.non.pair(1, 2, 1, 2)) + 
	(g22 * err.non.pair(2, 2, 1, 2))

s22 = 
	(g00 * err.non.pair(0, 0, 2, 2)) + 
	(g01 * err.non.pair(0, 1, 2, 2)) + 
	(g02 * err.non.pair(0, 2, 2, 2)) + 
	(g11 * err.non.pair(1, 1, 2, 2)) + 
	(g12 * err.non.pair(1, 2, 2, 2)) + 
	(g22 * err.non.pair(2, 2, 2, 2))



p00 = g00 * (
	( g00 * err.non.pair(0, 0, 0, 0)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 0, 0)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 0, 0)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 0, 0)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 0, 0)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 0, 0)  /  s22 ) )

p01 = g01 * (
	( g00 * err.non.pair(0, 0, 0, 1)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 0, 1)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 0, 1)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 0, 1)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 0, 1)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 0, 1)  /  s22 ) )

p02 = g02 * (
	( g00 * err.non.pair(0, 0, 0, 2)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 0, 2)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 0, 2)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 0, 2)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 0, 2)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 0, 2)  /  s22 ) )

p11 = g11 * (
	( g00 * err.non.pair(0, 0, 1, 1)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 1, 1)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 1, 1)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 1, 1)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 1, 1)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 1, 1)  /  s22 ) )

p12 = g12 * (
	( g00 * err.non.pair(0, 0, 1, 2)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 1, 2)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 1, 2)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 1, 2)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 1, 2)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 1, 2)  /  s22 ) )

p22 = g22 * (
	( g00 * err.non.pair(0, 0, 2, 2)  /  s00 ) + 
		( g01 * err.non.pair(0, 1, 2, 2)  /  s01 ) + 
		( g02 * err.non.pair(0, 2, 2, 2)  /  s02 ) + 
		( g11 * err.non.pair(1, 1, 2, 2)  /  s11 ) + 
		( g12 * err.non.pair(1, 2, 2, 2)  /  s12 ) + 
		( g22 * err.non.pair(2, 2, 2, 2)  /  s22 ) )


range(p00 + p01 + p02 + p11 + p12 + p22)




d = rbind(data.table(f = q, p = p00, e = g00, gt="00"),
					data.table(f = q, p = p01, e = g01, gt="01"),
					data.table(f = q, p = p02, e = g02, gt="02"),
					data.table(f = q, p = p11, e = g11, gt="11"),
					data.table(f = q, p = p12, e = g12, gt="12"),
					data.table(f = q, p = p22, e = g22, gt="22"))


ggplot(d) +  geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + geom_line(aes(x=f, y=e, colour = gt))







g00 = gen.ibd.pair(0, 0)
g01 = gen.ibd.pair(0, 1)
g02 = gen.ibd.pair(0, 2)
g11 = gen.ibd.pair(1, 1)
g12 = gen.ibd.pair(1, 2)
g22 = gen.ibd.pair(2, 2)


s00 = 
	(g00 * err.ibd.pair(0, 0, 0, 0)) + 
	(g01 * err.ibd.pair(0, 1, 0, 0)) + 
	(g02 * err.ibd.pair(0, 2, 0, 0)) + 
	(g11 * err.ibd.pair(1, 1, 0, 0)) + 
	(g12 * err.ibd.pair(1, 2, 0, 0)) + 
	(g22 * err.ibd.pair(2, 2, 0, 0))

s01 = 
	(g00 * err.ibd.pair(0, 0, 0, 1)) + 
	(g01 * err.ibd.pair(0, 1, 0, 1)) + 
	(g02 * err.ibd.pair(0, 2, 0, 1)) + 
	(g11 * err.ibd.pair(1, 1, 0, 1)) + 
	(g12 * err.ibd.pair(1, 2, 0, 1)) + 
	(g22 * err.ibd.pair(2, 2, 0, 1))

s02 = 
	(g00 * err.ibd.pair(0, 0, 0, 2)) + 
	(g01 * err.ibd.pair(0, 1, 0, 2)) + 
	(g02 * err.ibd.pair(0, 2, 0, 2)) + 
	(g11 * err.ibd.pair(1, 1, 0, 2)) + 
	(g12 * err.ibd.pair(1, 2, 0, 2)) + 
	(g22 * err.ibd.pair(2, 2, 0, 2))

s11 = 
	(g00 * err.ibd.pair(0, 0, 1, 1)) + 
	(g01 * err.ibd.pair(0, 1, 1, 1)) + 
	(g02 * err.ibd.pair(0, 2, 1, 1)) + 
	(g11 * err.ibd.pair(1, 1, 1, 1)) + 
	(g12 * err.ibd.pair(1, 2, 1, 1)) + 
	(g22 * err.ibd.pair(2, 2, 1, 1))

s12 = 
	(g00 * err.ibd.pair(0, 0, 1, 2)) + 
	(g01 * err.ibd.pair(0, 1, 1, 2)) + 
	(g02 * err.ibd.pair(0, 2, 1, 2)) + 
	(g11 * err.ibd.pair(1, 1, 1, 2)) + 
	(g12 * err.ibd.pair(1, 2, 1, 2)) + 
	(g22 * err.ibd.pair(2, 2, 1, 2))

s22 = 
	(g00 * err.ibd.pair(0, 0, 2, 2)) + 
	(g01 * err.ibd.pair(0, 1, 2, 2)) + 
	(g02 * err.ibd.pair(0, 2, 2, 2)) + 
	(g11 * err.ibd.pair(1, 1, 2, 2)) + 
	(g12 * err.ibd.pair(1, 2, 2, 2)) + 
	(g22 * err.ibd.pair(2, 2, 2, 2))



p00 = g00 * ( ( g00 * err.ibd.pair(0, 0, 0, 0)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 0, 0)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 0, 0)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 0, 0)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 0, 0)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 0, 0)  /  s22 ) )

p01 = g01 * ( ( g00 * err.ibd.pair(0, 0, 0, 1)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 0, 1)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 0, 1)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 0, 1)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 0, 1)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 0, 1)  /  s22 ) )

p02 = g02 * ( ( g00 * err.ibd.pair(0, 0, 0, 2)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 0, 2)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 0, 2)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 0, 2)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 0, 2)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 0, 2)  /  s22 ) )

p11 = g11 * ( ( g00 * err.ibd.pair(0, 0, 1, 1)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 1, 1)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 1, 1)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 1, 1)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 1, 1)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 1, 1)  /  s22 ) )

p12 = g12 * ( ( g00 * err.ibd.pair(0, 0, 1, 2)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 1, 2)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 1, 2)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 1, 2)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 1, 2)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 1, 2)  /  s22 ) )

p22 = g22 * ( ( g00 * err.ibd.pair(0, 0, 2, 2)  /  s00 ) + 
								( g01 * err.ibd.pair(0, 1, 2, 2)  /  s01 ) + 
								#( g02 * err.ibd.pair(0, 2, 2, 2)  /  s02 ) + 
								( g11 * err.ibd.pair(1, 1, 2, 2)  /  s11 ) + 
								( g12 * err.ibd.pair(1, 2, 2, 2)  /  s12 ) + 
								( g22 * err.ibd.pair(2, 2, 2, 2)  /  s22 ) )


range(p00 + p01 + p02 + p11 + p12 + p22, na.rm = T)




d = rbind(data.table(f = q, p = p00, e = g00, gt="00"),
					data.table(f = q, p = p01, e = g01, gt="01"),
					data.table(f = q, p = p02, e = g02, gt="02"),
					data.table(f = q, p = p11, e = g11, gt="11"),
					data.table(f = q, p = p12, e = g12, gt="12"),
					data.table(f = q, p = p22, e = g22, gt="22"))


ggplot(d) +  geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + geom_line(aes(x=f, y=e, colour = gt))




















p = (0:5000)/5000
q = 1-p

w = 0.01
r = 1 - w

e00 = r^2
e01 = 2 * w * r
e02 = w^2

e10 = 2 * w * r
e11 = 2 * ((w^2) + (r^2))
e12 = 2 * w * r

e20 = w^2
e21 = 2 * w * r
e22 = r^2


g0 = p^2
g1 = 2 * p * q
g2 = q^2


x0 = (g0 * e00) + (g1 * e10) + (g2 * e20)
x1 = (g0 * e01) + (g1 * e11) + (g2 * e21)
x2 = (g0 * e02) + (g1 * e12) + (g2 * e22)


y0 = g0 * (
	  ((g0 * e00) / x0) + 
		((g1 * e10) / x1) + 
		((g2 * e20) / x2) )

y1 = g1 * (
  	((g0 * e01) / x0) + 
		((g1 * e11) / x1) + 
		((g2 * e21) / x2) )

y2 = g2 * (
  	((g0 * e02) / x0) + 
		((g1 * e12) / x1) + 
		((g2 * e22) / x2) )

range(y0 + y1 + y2)

d = rbind(data.table(f = q, p = y0, e = g0, gt="0"),
					data.table(f = q, p = y1, e = g1, gt="1"),
					data.table(f = q, p = y2, e = g2, gt="2"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype proportion")



d = rbind(data.table(f = q, p = y0^2,        e = g0^2,         gt="00"),
					data.table(f = q, p = 2 * y0 * y1, e = 2 * g0 * g1,  gt="01"),
					data.table(f = q, p = 2 * y0 * y2, e = 2 * g0 * g2,  gt="02"),
					data.table(f = q, p = y1^2,        e = g1^2,         gt="11"),
					data.table(f = q, p = 2 * y1 * y2, e = 2 * g1 * g2,  gt="12"),
					data.table(f = q, p = y2^2,        e = g2^2,         gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")






p = (0:5000)/5000
q = 1-p

w = 0.01
r = 1 - w

e00 = r^2
e01 = 2 * w * r
e02 = w^2

e10 = 2 * w * r
e11 = 2 * ((w^2) + (r^2))
e12 = 2 * w * r

e20 = w^2
e21 = 2 * w * r
e22 = r^2


g0 = p^2
g1 = 2 * p * q
g2 = q^2

range(g0 + g1 + g2)


obs.g0 = (g0 * e00) + (g1 * e10/2) + (g2 * e20)
obs.g1 = (g0 * e01) + (g1 * e11/2) + (g2 * e21)
obs.g2 = (g0 * e02) + (g1 * e12/2) + (g2 * e22)

obs.p = sqrt(obs.g0)
obs.q = sqrt(obs.g2)


d = rbind(data.table(gt="00",  f = q,  e = g0,  p = obs.g0  ),
					data.table(gt="01",  f = q,  e = g1,  p = obs.g1  ),
					data.table(gt="02",  f = q,  e = g2,  p = obs.g2  ))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt))


x0 = (g0 * e00) + (g1 * e10) + (g2 * e20)
x1 = (g0 * e01) + (g1 * e11) + (g2 * e21)
x2 = (g0 * e02) + (g1 * e12) + (g2 * e22)



d = rbind(data.table(gt="00",  f = q,  e = g0,  p = x0  ),
					data.table(gt="01",  f = q,  e = g1,  p = x1  ),
					data.table(gt="02",  f = q,  e = g2,  p = x2  ))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt))



p00 = (g0 * e00) / x0
p10 = (g1 * e10) / x1
p20 = (g2 * e20) / x2

p01 = (g0 * e01) / x0
p11 = (g1 * e11) / x1
p21 = (g2 * e21) / x2

p02 = (g0 * e02) / x0
p12 = (g1 * e12) / x1
p22 = (g2 * e22) / x2

p0 = g0 * (p00 + p10 + p20)
p1 = g1 * (p01 + p11 + p21)
p2 = g2 * (p02 + p12 + p22)


obs.frq = sqrt(x2)
obs.p = sqrt(x0)
obs.q = sqrt(x2)

d = rbind(data.table(gt="00",  f = q,  e = p^4,            p = p0^2  ),
					data.table(gt="01",  f = q,  e = 4 * p^3 * q,    p = 2 * p0 * p1  ),
					data.table(gt="02",  f = q,  e = 2 * p^2 * q^2,  p = 2 * p0 * p2  ),
					data.table(gt="11",  f = q,  e = 4 * p^2 * q^2,  p = p1^2  ),
					data.table(gt="12",  f = q,  e = 4 * p * q^3,    p = 2 * p1 * p2  ),
					data.table(gt="22",  f = q,  e = q^4,            p = p2^2  ))

d = rbind(data.table(gt="00",  f = q,  e = p^3,                    p = p^3 + (p0^2 - g0^2)  ),
					data.table(gt="01",  f = q,  e = 2 * p^2 * q,            p = (2 * p^2 * q) + (2 * p0 * p1) - (2 * g0 * g1)  ),
					data.table(gt="02",  f = q,  e = 0,                      p = 0 + (2 * p0 * p2) - (2 * g0 * g2)  ),
					data.table(gt="11",  f = q,  e = (p^2 * q) + (p * q^2),  p = (p^2 * q) + (p * q^2) + p1^2 - g1^2  ),
					data.table(gt="12",  f = q,  e = 2 * p * q^2,            p = (2 * p * q^2) + (2 * p1 * p2) - (2 * g1 * g2)  ),
					data.table(gt="22",  f = q,  e = q^3,                    p = q^3 + (p2^2 - g2^2)  ))


spl = split(d, d$gt)
rng = rep(0, length(p))
for (gt in names(spl)) { rng = rng + spl[[gt]]$p }
range(rng)

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt))



x0 = (g0 * e00)^2  +  (2 * (g0 * e00) * (g1 * e10))  +  (2 * (g0 * e00) * (g2 * e20))  +  (g1 * e10)^2  +  (2 * (g1 * e10) * (g2 * e20))  +  (g2 * e20)^2
x1 = (g0 * e01)^2  +  (2 * (g0 * e01) * (g1 * e11))  +  (2 * (g0 * e01) * (g2 * e21))  +  (g1 * e11)^2  +  (2 * (g1 * e11) * (g2 * e21))  +  (g2 * e21)^2
x2 = (g0 * e02)^2  +  (2 * (g0 * e02) * (g1 * e12))  +  (2 * (g0 * e02) * (g2 * e22))  +  (g1 * e12)^2  +  (2 * (g1 * e12) * (g2 * e22))  +  (g2 * e22)^2

x0 = (g0 * e00) + (g1 * e10) + (g2 * e20)
x1 = (g0 * e01) + (g1 * e11) + (g2 * e21)
x2 = (g0 * e02) + (g1 * e12) + (g2 * e22)


p.00.00 = (g0 * e00) * (g0 * e00)     / (x0 * x0)
p.01.00 = 2 * (g0 * e00) * (g1 * e10) / (x0 * x1)
p.02.00 = 2 * (g0 * e00) * (g2 * e20) / (x0 * x2)
p.11.00 = (g1 * e10) * (g1 * e10)     / (x1 * x1)
p.12.00 = 2 * (g1 * e10) * (g2 * e20) / (x1 * x2)
p.22.00 = (g2 * e20) * (g2 * e20)     / (x2 * x2)

p.00.01 = (g0 * e00) * (g0 * e01)       / (x0 * x0)
p.01.01 = 2 * (g0 * e00) * (g1 * e11)   / (2 * x0 * x1)
p.02.01 = 2 * (g0 * e00) * (g2 * e21)   / (2 * x0 * x2)
p.11.01 = (g1 * e10) * (g1 * e11)       / (x1 * x1)
p.12.01 = 2 * (g1 * e10) * (g2 * e21)   / (2 * x1 * x2)
p.22.01 = (g2 * e20) * (g2 * e21)       / (x2 * x2)

p.00.02 = (g0 * e00) * (g0 * e02)       / (x0 * x0)
p.01.02 = 2 * (g0 * e00) * (g1 * e12)   / (2 * x0 * x1)
p.02.02 = 2 * (g0 * e00) * (g2 * e22)   / (2 * x0 * x2)
p.11.02 = (g1 * e10) * (g1 * e12)       / (x1 * x1)
p.12.02 = 2 * (g1 * e10) * (g2 * e22)   / (2 * x1 * x2)
p.22.02 = (g2 * e20) * (g2 * e22)       / (x2 * x2)


d = rbind(data.table(gt="00",  f = q,  e = p^4,            p = p^4 * (p.00.00 + p.01.00 + p.02.00 + p.11.00 + p.12.00 + p.22.00)  ),
					data.table(gt="01",  f = q,  e = 4 * p^3 * q,    p = 4 * p^3 * q * (p.00.01 + p.01.01 + p.02.01 + p.11.01 + p.12.01 + p.22.01)  ),
					data.table(gt="02",  f = q,  e = 2 * p^2 * q^2,  p = 2 * p^2 * q^2 * (p.00.02 + p.01.02 + p.02.02 + p.11.02 + p.12.02 + p.22.02)  ),
					data.table(gt="11",  f = q,  e = 4 * p^2 * q^2,  p = (p01 + p11 + p21)^2  ),
					data.table(gt="12",  f = q,  e = 4 * p * q^3,    p = 2 * (p01 + p11 + p21) * (p02 + p12 + p22)  ),
					data.table(gt="22",  f = q,  e = q^4,            p = (p02 + p12 + p22)^2  ))

spl = split(d, d$gt)
rng = rep(0, length(p))
for (gt in names(spl)) { rng = rng + spl[[gt]]$p }
range(rng)

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22", size = 1) + 
	geom_line(aes(x=f, y=e, colour = gt))




d = rbind(data.table(gt="00",  f = q,  e = p^3,                    p = (p00 + p10 + p20) * sqrt(p00 + p10 + p20)  ),
					data.table(gt="01",  f = q,  e = 2 * p^2 * q,            p = (p01 + p11 + p21) * sqrt(p00 + p10 + p20)  ),
					data.table(gt="02",  f = q,  e = 0,                      p = ((p00 + p10 + p20) * (p02 + p12)) + (p10 + p20) * (p02 + p12 + p22) ),
					data.table(gt="11",  f = q,  e = (p^2 * q) + (p * q^2),  p = ((p01 + p11 + p21) * sqrt(p00 + p10 + p20) / 2) + ((p01 + p11 + p21) * sqrt(p02 + p12 + p22) / 2) ),
					data.table(gt="12",  f = q,  e = 2 * p * q^2,            p = (p01 + p11 + p21) * sqrt(p02 + p12 + p22)  ),
					data.table(gt="22",  f = q,  e = q^3,                    p = (p02 + p12 + p22) * sqrt(p02 + p12 + p22)  ))

spl = split(d, d$gt)
rng = rep(0, length(p))
for (gt in names(spl)) { rng = rng + spl[[gt]]$p }
range(rng)

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt))


y0 = g0 * (
	((g0 * e00) / x0) + 
  ((g1 * e10) / x1) + 
  ((g2 * e20) / x2) )

y1 = g1 * (
	((g0 * e01) / x0) + 
  ((g1 * e11) / x1) + 
  ((g2 * e21) / x2) )

y2 = g2 * (
	((g0 * e02) / x0) + 
  ((g1 * e12) / x1) + 
  ((g2 * e22) / x2) )

range(y0 + y1 + y2)

d = rbind(data.table(f = q, p = y0, e = g0, gt="0"),
					data.table(f = q, p = y1, e = g1, gt="1"),
					data.table(f = q, p = y2, e = g2, gt="2"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype proportion")



d = rbind(data.table(f = q, p = y0 * sqrt(y0),        e = p^3,          gt="00"),
					data.table(f = q, p = y1 * sqrt(y0), e = 2 * p^2 * q,  gt="01"),
					data.table(f = q, p = y2 * y0, e = 0,            gt="02"),
					data.table(f = q, p = y1^2,        e = (p^2 * q) + (p * q^2), gt="11"),
					data.table(f = q, p = 2 * y1 * y2, e = 2 * p * q^2,  gt="12"),
					data.table(f = q, p = y2 * sqrt(y2),        e = q^3,          gt="22"))

ggplot(d) +  
	geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + 
	geom_line(aes(x=f, y=e, colour = gt)) +
	theme(aspect.ratio = 1,
				legend.title = element_blank()) +
	xlab("Allele frequency") + 
	ylab("Genotype pair proportion")






gt0 = p^4 + (2 * p^3 * q) + (1 * p^2 * q^2)
gt1 = (4 * p^2 * q^2) + (2 * p^3 * q) + (2 * q^3 * p)
gt2 = q^4 + (2 * q^3 * p) + (1 * q^2 * p^2)

range(gt0 + gt1 + gt2)

d = rbind(data.table(f = q, p = gt0, gt = "0"),
					data.table(f = q, p = gt1, gt = "1"),
					data.table(f = q, p = gt2, gt = "2"))

ggplot(d) + geom_line(aes(x = f, y = p, colour = gt))




gt0 = p^3 + (p^2 * q)
gt1 = (p^2 * q) + (p * q^2) + (p^2 * q) + (q^2 * p)
gt2 = q^3 + (q^2 * p)

range(gt0 + gt1 + gt2)

d = rbind(data.table(f = q, p = gt0, gt = "0"),
					data.table(f = q, p = gt1, gt = "1"),
					data.table(f = q, p = gt2, gt = "2"))

ggplot(d) + geom_line(aes(x = f, y = p, colour = gt))














ep = sqrt(y0)
eq = sqrt(y2)

d = rbind(data.table(f = q, p = y0^1.5,   e = g0^1.5,      gt="00"),
					data.table(f = q, p = 2 * y0 * y2^0.5,              e = 2 * p^2 * q,            gt="01"),
					data.table(f = q, p = 0,                        e = 0,                      gt="02"),
					data.table(f = q, p = (ep^2 * eq) + (ep * eq^2),  e = (p^2 * q) + (p * q^2),  gt="11"),
					data.table(f = q, p = 2 * ep * eq^2,           e = 2 * p * q^2,            gt="12"),
					data.table(f = q, p = eq^3,                     e = q^3,                    gt="22"))

ggplot(d) +  geom_line(aes(x=f, y=p, colour = gt), linetype = "22") + geom_line(aes(x=f, y=e, colour = gt))





d = rbind(data.table(f = q, p = (g0 * e00) / ((g0 * e00) + (g1 * e01) + (g2 * e02)), from="0", to = "0"),
					data.table(f = q, p = (g0 * e10) / ((g0 * e10) + (g1 * e11) + (g2 * e12)), from="1", to = "0"),
					data.table(f = q, p = (g0 * e20) / ((g0 * e20) + (g1 * e21) + (g2 * e22)), from="2", to = "0"),
					
					data.table(f = q, p = (g1 * e01) / ((g0 * e00) + (g1 * e01) + (g2 * e02)), from="0", to = "1"),
					data.table(f = q, p = (g1 * e11) / ((g0 * e10) + (g1 * e11) + (g2 * e12)), from="1", to = "1"),
					data.table(f = q, p = (g1 * e21) / ((g0 * e20) + (g1 * e21) + (g2 * e22)), from="2", to = "1"),
					
					data.table(f = q, p = (g2 * e02) / ((g0 * e00) + (g1 * e01) + (g2 * e02)), from="0", to = "2"),
					data.table(f = q, p = (g2 * e12) / ((g0 * e10) + (g1 * e11) + (g2 * e12)), from="1", to = "2"),
					data.table(f = q, p = (g2 * e22) / ((g0 * e20) + (g1 * e21) + (g2 * e22)), from="2", to = "2"))

ggplot(d) + facet_grid(to~from) + geom_line(aes(x=f, y=p, colour = from))






f = function(i, j) {
	
	g = choose(2, 2 - (i + j))
	
	a = choose(2, i)
	b = choose(2, j)
	
	g*a*b
}















