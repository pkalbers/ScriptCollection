

source("MUTage_functions.R")

load("rea2.RData")



#Get positions of mutations
pos<-real_dat[, 1]


#Set Ne
N0 <- 10000 



#library(gtools)
#library(Matrix)

#setwd( "U:/Maurice");

do.setup<-FALSE;

if (do.setup) {
	cat("\nReading data");
	load("rea2.RData");
	cat("\nDone");

	#Get positions of mutations
	pos<-real_dat[,1];

	#Read in functions
	source("MUTage_functions.R")

	#Set Ne
	N0 <- 10000 

	#Parameters for file handling
	block.size<-10000;

	#Reorganise input files
	d<-real_dat[,-1];
	rm(real_dat);

	#Make genetic map
	gmap<-approx(tav[,1],tav[,3], xout=pos, rule=2)$y*4*N0/100;

	#Decide whether to include uncertainty about rec breaks
	hard.r.bks<-FALSE;

	
	#Data params
	l<-nrow(d);
	n<-ncol(d);
	k<-100;  ##k<-40
	pairs.max<-10000;
	alpha_n<-sum(1/(1:n));
	win_max<-2000;

	#Grid for times (in 2Ne units)
	tt.seq<-seq(0,5,l=1001);
	tt.seq[1]<-tt.seq[2]/2;
	tt.pos<-rep(0, length(tt.seq));

	f<-rowSums(d);

	#Stuff for breakpoint function
	tmp<-f/n;
	ff<-tmp^2+(1-tmp)^2;
	ff<-cumsum(log(ff));
	hh<-log(2*tmp*(1-tmp));

	#Mutation rate and N_e
	mut.rate.per.gen<-1.2e-8;
	n_e<-1e4;
	theta.per.site<-4*n_e*mut.rate.per.gen;         #Ne is just used as an internal scaling factor
	theta_est_Watterson<-l/(alpha_n*diff(range(pos)));
	theta_est_pwd<-sum(2*f/n*(1-f/n))/diff(range(pos));

	#Some storage space
	op<-array(0, c(pairs.max, 10));
	colnames(op)<-c("ID1", "ID2", "Mut.1", "Mut.2", "S.L", "S.R", "BK.L", "BK.R", "last.S.L", "last.S.R");

}




#Choose sites to analyse
pp<-50001:50002;


pt.estimate<-array(0, c(length(pp), 8));
colnames(pt.estimate)<-c("Post.mean", "Post.mode", "Post.median", "Post.2.5", "Post.97.5", "Robust", "Site", "Freq");
pt.estimate[,7]<-pos[pp];
pt.estimate[,8]<-f[pp];

ct.var<-0;

#Now iterate over sites
for (p in pp) if (f[p]>0 & f[p]<n) {
	
	cat("\nDoing site ", p-pp[1]+1, " of ", length(pp), sep="");
	
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
	#   tmp3[,2]<-sample(wni, l.i, replace=TRUE);           #Currently just a random selection
	for (i in 1:l.i) {
		if (p<l) {
			dp<-apply(d[(p+1):win.high,wni,drop=F]!=d[(p+1):win.high,tmp3[i,1]], 2, find.first.true.element);
		} else {dp<-rep(0, length(wni));}
		if (p>1) {
			dm<-apply(d[(p-1):win.low,wni,drop=F]!=d[(p-1):win.low,tmp3[i,1]], 2, find.first.true.element);
		} else {dm<-rep(0, length(wni));}
		tmp3[i,2]<-wni[which.max(dp+dm)];
	}
	
	#Count pairs
	n.k<-nrow(tmp1)+nrow(tmp2)+nrow(tmp3);
	ii<-1:n.k;
	
	
	
	#Start filling op
	op[ii,1]<-c(tmp1[,1], tmp2[,1], tmp3[,1]);
	op[ii,2]<-c(tmp1[,2], tmp2[,2], tmp3[,2]);
	
	#Randomise
	#for (i in 1:n.k) {
	#op[i,1:2]<-sample(n, 2, replace=F);
	#}
	
	op[ii,3]<-as.vector(d[p,op[ii,1]]);
	op[ii,4]<-as.vector(d[p,op[ii,2]]);
	op[ii,5:10]<-0;
	
	
	
	
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
		last.s<-p;
		while(flag) {
			j<-j+1;
			f.within<-sum(d[b.r[j],wi]);
			if ((f.within < f.focal) & (f.within < f[b.r[j]])) {
				flag<-FALSE;
				break.r<-b.r[j];
			} else {ct<-ct+1; last.s<-b.r[j];}
			if (j==length(b.r) & flag) {flag<-FALSE; break.r<-l+1;}
		}
		op[i,8]<-break.r;	#This is the location of the first difference - i.e. the point at which the break is detected
		op[i,6]<-ct;
		op[i,10]<-last.s;
		if (op[i,4]==0) op[i,6]<-op[i,6]+1;	#Add +1 for focal SNP if discordant
		
		#Analyse to LHS of pos
		b.l<-tmp[tmp<p];
		flag<-TRUE;
		if (length(b.l)==0) {flag<-FALSE; break.l<-0;}
		j<-length(b.l)+1;
		ct<-0;
		last.s<-p;
		while(flag) {
			j<-j-1;
			f.within<-sum(d[b.l[j],wi]);
			if ((f.within < f.focal) & (f.within < f[b.l[j]])) {
				flag<-FALSE;
				break.l<-b.l[j];
			} else {ct<-ct+1; last.s<-b.l[j];}
			if (j==1 & flag) {flag<-FALSE; break.l<-0;}
		}
		op[i,7]<-break.l;	#This is the location of the first difference - i.e. the point at which the break is detected
		op[i,5]<-ct;
		op[i,9]<-last.s;
	}
	
	
	#Now infer age distribution
	
	if (hard.r.bks) {
		a.params<-1+op[ii,5]+op[ii,6]+as.integer(op[ii,7]>0)+as.integer(op[ii,8]<=l);
		#    a.params<-1+op[ii,5]+op[ii,6];
		tmp1<-op[ii,8];
		tmp1[tmp1>l]<-l;
		tmp2<-op[ii,7];
		tmp2[tmp2<1]<-1;
		mut.length<-(pos[tmp1]-pos[tmp2]+1)*theta.per.site;
		rec.length<-(gmap[tmp1]-gmap[tmp2]);
		b.params<-1+mut.length+rec.length;
		#    b.params<-1+mut.length;
	}
	
	
	tt.pos[]<-0;
	plot(tt.seq,rep(0, length(tt.seq)),type="n", ylim=c(0,1), xlab="Time", ylab="CCF", log="x", main=paste("Pos = ", pos[p], ": f = ", f[p], sep=""));
	v.low<-0;
	v.high<-0;
	for (i in 1:n.k) {
		
		op.i<-op[i,];
		position<-p;
		
		if (!hard.r.bks) {
			llk.surf.i<-llk.tmrca.pair(op[i,], position=p, t.grid=tt.seq, gmap=gmap, pos=pos, ff=ff, hh=hh, 
																 theta.per.site=theta.per.site, debug.flag=FALSE);
			mx<-max(llk.surf.i);
			tmp1<-exp(llk.surf.i-mx);
			tmp1<-tmp1/sum(tmp1);
			cdf.i<-cumsum(tmp1);
			cdf.i<-1e-8+(1-2e-8)*cdf.i
		}
		
		if (op[i,4]==1) {
			if (hard.r.bks) {
				tt.pos<-tt.pos+pgamma(tt.seq, a.params[i], b.params[i], log=T, lower=T);
				lines(tt.seq, pgamma(tt.seq, a.params[i], b.params[i], lower=T), col="blue");
			} else {
				tt.pos<-tt.pos+log(cdf.i);
				lines(tt.seq, cdf.i, col="blue");
				v.low<-v.low+log(tt.seq[which.min(abs(cdf.i-0.5))]);
			}
		}
		if (op[i,4]==0) {
			if (hard.r.bks) {
				tt.pos<-tt.pos+pgamma(tt.seq, a.params[i], b.params[i], log=T, lower=F);
				lines(tt.seq, pgamma(tt.seq, a.params[i], b.params[i], lower=F), col="red");
			} else {
				tt.pos<-tt.pos+log(1-cdf.i);
				lines(tt.seq, 1-cdf.i, col="red");
				v.high<-v.high+log(tt.seq[which.min(abs(cdf.i-0.5))]);
			}
		}
		if (sum(is.na(tt.pos))>0) break();
	}
	mx<-max(tt.pos);
	tt.pos<-exp(tt.pos-mx);
	tt.pos<-tt.pos/sum(tt.pos);
	#plot(tt.seq, tt.pos, log="x", xlab="Age", ylab="Relative Posterior", type="l", 
	#    col="blue", main=paste("Pos = ", p, ": f = ", f[p], sep=""));
	lines(tt.seq, tt.pos/max(tt.pos), type="l", col="black", lwd=2);
	#abline(v=true.age[p], lty="dotted", col="red");
	
	n.within<-sum(op[ii,4]);
	if (n.within>0) v.low<-v.low/n.within;
	v.high<-v.high/(max(n.k-n.within, 1));
	
	if (n.within>0) {con.est<-(v.low+v.high)/2;} else {con.est<-v.high - log(2);}
	abline(v=exp(c(v.low, v.high, con.est)), col=c("blue", "red", "lightblue"));
	pt.estimate[ct.var, 6]<-con.est;
	
	
	#Get pseudo-95% ETPI
	tt.cum<-cumsum(tt.pos);
	oo<-approx(x=c(0,tt.cum), y=c(0,tt.seq), xout=c(0.025, 0.975));
	
	pt.estimate[ct.var,1]<-sum(tt.seq*tt.pos);
	pt.estimate[ct.var,2]<-tt.seq[which.max(tt.pos)];
	pt.estimate[ct.var,3]<-tt.seq[which.min(abs(cumsum(tt.pos)-0.5))];
	
	pt.estimate[ct.var,4:5]<-oo$y;  
	
}




write.table(pt.estimate, file="output.txt", quote=F, row=F, col=T);



