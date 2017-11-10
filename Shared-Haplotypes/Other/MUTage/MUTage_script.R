
#Code to estimate mutation age for a few thousand samples

is.desktop<-FALSE;

#Setwd and read in functions
if (is.desktop) setwd("U:/Maurice");
if (!is.desktop) setwd("C:/Users/mcvean/Documents/Maurice");
source("MUTage_functions.R");

sample.safe<-function(ids, k, replace=FALSE) {
	l.i<-length(ids);
	if (l.i==1) {
		return(ids);
	} else {
		return(sample(ids, min(k, l.i), replace=replace));
	}
}

find.first.true.element<-function(vtr) {
	tmp<-which(vtr);
	if (length(tmp)==0) { 
		return(length(vtr)+1);
	} else {return(tmp[1])};
}

#Parameters for file handling
block.size<-10000;

#Load data
if (!("d" %in% ls())) {
	load("kde-s.R");
	d<-kde[,1:(ncol(kde)-2)];				#Haplotypes only - as 0/1
	pos<-kde[,ncol(kde)-1];					#Position of each SNP
	age<-kde[,ncol(kde)];					#Age of each SNP in generations
	gmap<-approx(oo[,1], oo[,3], xout=pos);		#Genetic map position (cM) for each SNP
	rm(kde);
	rm(oo);
}

#Data params
l<-nrow(d);
n<-ncol(d);
k<-40;
pairs.max<-1000;
alpha_n<-sum(1/(1:n));
win_max<-1000;

#Grid for times (in 2Ne units)
tt.seq<-seq(0,5,l=1001);
tt.seq[1]<-tt.seq[2]/2;
tt.pos<-rep(0, length(tt.seq));

#Mutation rate and N_e
mut.rate.per.gen<-1e-8;
n_e<-1e4;
theta.per.site<-4*n_e*mut.rate.per.gen;			#Ne is just used as an internal scaling factor
theta_est_Watterson<-l/(alpha_n*diff(range(pos)));
theta_est_pwd<-sum(2*f/n*(1-f/n))/diff(range(pos));

#Some storage space
op<-array(0, c(pairs.max, 7));
colnames(op)<-c("ID1", "ID2", "Mut.1", "Mut.2", "S", "OR.L", "OR.R");


#Calculate frequencies
f<-rep(0,l);
flag<-TRUE;
ct<-0;
while(flag) {
	cat("\rCurrently start location = ", ct+1, "\t\t");
	i.low<-(ct+1);
	i.high<-min(ct+block.size, l);
	f[i.low:i.high]<-apply(d[i.low:i.high,], 1, sum);
	ct<-ct+block.size;
	if (ct>=l) flag<-FALSE;
}



#Choose sites to analyse
pp<-c(10000:10200);
pt.estimate<-array(0, c(length(pp), 3));
colnames(pt.estimate)<-c("Post.mean", "Post.mode", "Post.median");
ct.var<-0;

#Now iterate over sites
for (p in pp) if (f[p]>0 & f[p]<n) {

	ct.var<-ct.var+1;
	win.low<-max(1, p-win_max);
	win.high<-min(l, p+win_max);

	#Work out who has mutation
	wi<-which(d[p,]==1);
	wni<-setdiff(1:n, wi);
	n.i<-length(wi);

	#Sample pairs of haplotypes to analyse
	#Within clade
	if (n.i>1) {
		tmp1<-sample.pairs(wi, k);
	} else tmp1<-array(0, c(0,2));

	#Across clade
	tmp2<-sample.pairs2(wi, wni, k);

	#Augment pairs with maximal sharing outside clade
	who.to.analyse<-sample.safe(wi, k);
	l.i<-length(who.to.analyse);
	tmp3<-array(0, c(l.i, 2));
	tmp3[,1]<-who.to.analyse;
#	tmp3[,2]<-sample(wni, l.i, replace=TRUE);			#Currently just a random selection
	for (i in 1:l.i) {
		dp<-apply(d[(p+1):win.high,wni]!=d[(p+1):win.high,tmp3[i,1]], 2, find.first.true.element);
		dm<-apply(d[(p-1):win.low,wni]!=d[(p-1):win.low,tmp3[i,1]], 2, find.first.true.element);
		tmp3[i,2]<-wni[which.max(dp+dm)];
	}
	

	#Count pairs
	n.k<-nrow(tmp1)+nrow(tmp2)+nrow(tmp3);
	ii<-1:n.k;
	
	#Start filling op
	op[ii,1]<-c(tmp1[,1], tmp2[,1], tmp3[,1]);
	op[ii,2]<-c(tmp1[,2], tmp2[,2], tmp3[,2]);
	op[ii,3]<-d[p,op[ii,1]];
	op[ii,4]<-d[p,op[ii,2]];
	op[ii,5:7]<-0;

	#Now go through each pair getting breaks and number of mutations
	for (i in 1:n.k) {

		tmp<-which(d[,op[i,1]] != d[,op[i,2]]);
		f.focal<-f[p];

		#Analyse to RHS of pos
		b.r<-tmp[tmp>p];
		flag<-TRUE;
		if (length(b.r)==0) {flag<-FALSE; break.r<-l+1;}
		j<-0;
		ct<-0;
		while(flag) {
			j<-j+1;
			f.within<-sum(d[b.r[j],wi]);
			if ((f.within < f.focal) & (f.within < f[b.r[j]])) {
				flag<-FALSE;
				break.r<-b.r[j];
			} else {ct<-ct+1;}
			if (j==length(b.r) & flag) {flag<-FALSE; break.r<-l+1;}
		}
		op[i,7]<-break.r;
		op[i,5]<-ct;

		#Analyse to LHS of pos
		b.l<-tmp[tmp<p];
		flag<-TRUE;
		if (length(b.l)==0) {flag<-FALSE; break.l<-0;}
		j<-length(b.l)+1;
		ct<-0;
		while(flag) {
			j<-j-1;
			f.within<-sum(d[b.l[j],wi]);
			if ((f.within < f.focal) & (f.within < f[b.l[j]])) {
				flag<-FALSE;
				break.l<-b.l[j];
			} else {ct<-ct+1;}
			if (j==1 & flag) {flag<-FALSE; break.l<-0;}
		}
		op[i,6]<-break.l;
		op[i,5]<-op[i,5]+ct;
	}

	#Now infer age distribution
	a.params<-1+op[ii,5]+as.integer(op[ii,6]>0)+as.integer(op[ii,7]<=l);
	tmp1<-op[ii,7];
	tmp1[tmp1>l]<-l;
	tmp2<-op[ii,6];
	tmp2[tmp2<1]<-1;
	mut.length<-(pos[tmp1]-pos[tmp2])*theta.per.site;
	rec.length<-0.04*n_e*(gmap$y[tmp1]-gmap$y[tmp1]);
	b.params<-1+mut.length+rec.length;

	tt.pos[]<-0;
	for (i in 1:n.k) {
		if (op[i,4]==1) tt.pos<-tt.pos+pgamma(tt.seq, a.params[i], b.params[i], log=T, lower=T);
		if (op[i,4]==0) tt.pos<-tt.pos+pgamma(tt.seq, a.params[i], b.params[i], log=T, lower=F);
	}
	mx<-max(tt.pos);
	tt.pos<-exp(tt.pos-mx);
	tt.pos<-tt.pos/sum(tt.pos);
	plot(tt.seq, tt.pos, log="x", xlab="Age", ylab="Relative Posterior", type="l", 
		col="blue", main=paste("Pos = ", p, ": f = ", f[p], sep=""));
	abline(v=age[p], lty="dotted", col="red");

	pt.estimate[ct.var,1]<-sum(tt.seq*tt.pos);
	pt.estimate[ct.var,2]<-tt.seq[which.max(tt.pos)];
	pt.estimate[ct.var,3]<-tt.seq[which.min(abs(cumsum(tt.pos)-0.5))];
	
}

plot(age[pp], pt.estimate[,1], pch=19, col="blue", xlab="True Age", ylab="Inferred Age", log="xy");
abline(0,1,col="red", lty="dotted");
cat("\nCor true age with estimated age = ", cor(age[pp], pt.estimate, method="spearman"));
cat("\nCor true age with freq = ", cor(age[pp], f[pp], method="spearman"));


