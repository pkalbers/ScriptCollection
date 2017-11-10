

#Functions for mutation age

######################################################################
#Function to sample overcoming inbuilt R issue when lengths(ids)==1
######################################################################

sample.safe<-function(ids, k, replace=FALSE) {
	l.i<-length(ids);
	if (l.i==1) {
		return(ids);
	} else {
		return(sample(ids, min(k, l.i), replace=replace));
	}
}

#########################################################
#Funciton to first first true element in a Boolean array
#########################################################

find.first.true.element<-function(vtr) {
	tmp<-which(vtr);
	if (length(tmp)==0) { 
		return(length(vtr)+1);
	} else {return(tmp[1])};
}


###############################################
#Function returns up to n non-redundant pairs 
###############################################

sample.pairs<-function(ids, n, debug=FALSE) {
	
	ni<-length(ids);
	nmax<-ni*(ni-1)/2;
	if (n>nmax) n<-nmax;
	
	op<-array(0, c(n, 2));
	wp<-sample(nmax, n);
	
	if (!debug) {
		
		ii<-1:(ni-1);
		ll<-(ii-1)*(2*ni-ii)/2+1;
		op[,1]<-findInterval(wp, ll);
		op[,2]<-wp-ll[op[,1]]+1+op[,1];
	} else {
		tmp<-array(0, c(nmax, 3));
		ct<-0;
		for (i in 1:(ni-1)) for (j in (i+1):ni) {
			ct<-ct+1;
			tmp[ct,]<-c(ct,i,j);
		}
		op<-tmp[wp,2:3];
	}
	
	op[,1]<-ids[op[,1]];
	op[,2]<-ids[op[,2]];
	return(op);
}

######################################################
#Function to sample pairs for list1xlist2
######################################################

sample.pairs2<-function(ids1, ids2, n, debug=FALSE) {
	
	ni<-length(ids1);
	nj<-length(ids2);
	nmax<-ni*nj;
	if (n>nmax) n<-nmax;
	
	op<-array(0, c(n,2));
	wp<-sample(nmax, n);
	
	if (!debug) {
		ii<-1:ni;
		ll<-(ii-1)*nj+1;
		op[,1]<-findInterval(wp, ll);
		op[,2]<-wp-ll[op[,1]]+1;
	} else {
		tmp<-array(0, c(nmax, 3));
		ct<-0;
		for (i in 1:ni) for (j in 1:nj) {
			ct<-ct+1;
			tmp[ct,]<-c(ct, i, j);
		}
		op<-tmp[wp,2:3];
	}
	
	op[,1]<-ids1[op[,1]];
	op[,2]<-ids2[op[,2]];
	
	return(op);
}


########################################################################################################
#Function to calculate empirical probability of breaking a haplotype at each position, left and right
#ebpd = empirical breakpoint distribution
########################################################################################################

ebpd<-function(seq, position=1, n.pairs=1000, naive="TRUE", do.plot="TRUE") {
	
	n<-ncol(seq);
	n.pairs<-min(n.pairs, n*(n-1)/2);
	pairs<-sample.pairs(1:n, n.pairs)
	scr<-array(0, c(n.pairs, 2));
	mx<-nrow(seq)+1;
	
	for (i in 1:n.pairs) {
		del<-which(d[,pairs[i,1]] != d[,pairs[i,2]]);
		scr[i,1]<-max(c(0,del[del<position]));
		scr[i,2]<-min(c(del[del>position], mx));
	}
	
	ecdf.l<-ecdf(-scr[,1]);
	ecdf.r<-ecdf(scr[,2]);
	
	if (do.plot) {
		par(mfrow=c(1,2));
		plot(ecdf.l, xlab="Position rel to focal", ylab="CDF", col="blue", main=paste("Position", position, ": Left", sep=" "));
		plot(ecdf.r, xlab="Position rel to focal", ylab="CDF", col="blue", main=paste("Position", position, ": Right", sep=" "));
		par(mfrow=c(1,1));
	}
	
	#Compare to naive function
	if (naive) {	
		
		#Left
		l.mn<-min(scr[,1]);
		f.win<-rowSums(seq[l.mn:position,]);
		v1<-ecdf.l(-((position-1):l.mn));
		hh<-rev(2*f.win[-length(f.win)]/n*(1-f.win[-length(f.win)]/n));
		v2<-hh*exp(cumsum(c(0,log(1-hh[-length(hh)]))));
		
		par(mfrow=c(1,2));
		plot(v1, type="l", xlab="Position rel to focal", log="x", ylim=c(0,1),
				 ylab="CDF", col="blue", main=paste("Position", position, ": Left", sep=" "));
		lines(cumsum(v2), col="red");
		
		#Right
		r.mx<-max(scr[,2]);
		f.win<-rowSums(seq[position:r.mx,]);
		v1<-ecdf.r((position+1):r.mx);
		hh<-2*f.win[-1]/n*(1-f.win[-1]/n);
		v2<-hh*exp(cumsum(c(0,log(1-hh[-length(hh)]))));
		
		plot(v1, type="l", xlab="Position rel to focal", log="x", ylim=c(0,1), 
				 ylab="CDF", col="blue", main=paste("Position", position, ": Right", sep=" "));
		lines(cumsum(v2), col="red");
		
		par(mfrow=c(1,1));
	}
	
	return(list(ecdf.l=ecdf.l, ecdf.r=ecdf.r));
	
}



###############################################################################
#Function to calculate freq-based approximation to break-point calculation
#ff is a cumulative map of log homozygosities
#hh is a per SNP list of log heterozygosities
#Note may have list of b for any segbreak
###############################################################################

approx.bpf<-function(span, segbreak, cum_log_hom, log_het) {
	
	l.sq<-length(cum_log_hom);
	
	if (segbreak %in% span) {
		cat("\n\nError in BPF - cannot have break and detected break in same position\n\n");
		return();
	}
	
	#Calculating prob break at segbreak and cope with boundary cases
	if (segbreak > 0 && segbreak <= l.sq) {
		v1 <- rep(log_het[segbreak], length(span));
	} else {
		v1 <- rep(0, length(span));
	}
	
	#Calculate sum of log identity up to point - coping with boundary cases
	if (segbreak > span[1]) {		#Right hand calculation
		
		if (segbreak <= l.sq) {
			v1 <- v1 + cum_log_hom[segbreak-1] - cum_log_hom[span]; # not boundary
		} else {
			v1 <- v1 + cum_log_hom[l.sq] - cum_log_hom[span]; # boundary!
		}
		
	} else {			#LH calculation
		
		if (segbreak > 0) {
			v1 <- v1 + cum_log_hom[span - 1] - cum_log_hom[segbreak]; # not boundary
		} else {
			wi <- which(span > 1); # exclude boundary
			v1[wi] <- v1[wi] + cum_log_hom[span[wi] - 1]; # boundary!
		}
		
	}
	
	return(v1);
}



#####################################################################
#Function to calculate llk surface for tmrca given S and breakpoints
#Sums over unknown location of breakpoint
#####################################################################

llk.tmrca.pair<-function(op.i, position, t.grid, gmap, pos, ff, hh, theta.per.site=0.001, debug.flag=FALSE, include.tij.prior=TRUE, max.poss.bks=1000000) {
	
	s.l<-op.i[5];
	s.r<-op.i[6];
	bk.l<-op.i[7];
	bk.r<-op.i[8];
	last.s.l<-op.i[9];
	last.s.r<-op.i[10];
	
	l<-length(pos);
	
	
	#RHS
	#First get P(bk.r | b) for b in {last.s.r, bk.r - 1);
	#Note break means in interval i to i+1	
	
	poss.bks.r <- max(position, bk.r - max.poss.bks):(bk.r-1);
	
	l.b.r<-length(poss.bks.r);
	#if (debug.flag) {poss.bks.r<-c(bk.r-1); l.b.r<-1;}
	p.bk.r <- approx.bpf(poss.bks.r, bk.r, ff, hh);
	
	#Work out physical lengths - taking edge cases into account
	if (poss.bks.r[l.b.r] >= l) {
		lengths.r <- rep(0, l.b.r);
		lengths.r[l.b.r] <- pos[l] - pos[position] + 1e-8;
		lengths.r[1:(l.b.r-1)] <- (pos[poss.bks.r[-l.b.r]] + pos[poss.bks.r[-l.b.r] + 1]) / 2 - pos[position] + 1e-8;
	} else {
		lengths.r <- (pos[poss.bks.r] + pos[poss.bks.r + 1]) / 2 - pos[position] + 1e-8;
	}
	
	#Work out genetic lengths of region unbroken by rec
	rec.lengths.r <- gmap[poss.bks.r] - gmap[position] + 1e-8;
	
	#Work out amount of recombination in the interval with break taking edge cases into account
	if (poss.bks.r[l.b.r]>=l) {
		del.rec.r<-rep(0, l.b.r);
		del.rec.r[1:(l.b.r-1)] <- gmap[poss.bks.r[-l.b.r]+1] - gmap[poss.bks.r[-l.b.r]] + 1e-8;
		del.rec.r[l.b.r] <- 1e-8;
	} else {
		del.rec.r <- gmap[poss.bks.r+1] - gmap[poss.bks.r] + 1e-8;
	}
	
	
	#LHS
	#Now get P(bk.l | b) for b in {bk.l + 1 : last.s.l);
	#Note that to LHS, the definition of break changes to break in i-1 to i interval
	
	poss.bks.l <- (bk.l + 1) : min(position, bk.l + max.poss.bks);
	
	l.b.l<-length(poss.bks.l);
	#if (debug.flag) {poss.bks.l<-c(bk.l+1); l.b.l<-1;}
	p.bk.l <- approx.bpf(poss.bks.l, bk.l, ff, hh);
	
	#Work out physical lengths - taking edge cases into account
	if (poss.bks.l[1] == 1) {
		lengths.l <- rep(0,l.b.l);
		lengths.l[1] <- pos[position] - pos[1] + 1e-8;
		lengths.l[2:l.b.l] <- pos[position] - (pos[poss.bks.l[-1]] + pos[poss.bks.l[-1]-1]) / 2 + 1e-8;
	} else {
		lengths.l <- pos[position] - (pos[poss.bks.l] + pos[poss.bks.l-1]) / 2 + 1e-8;
	}
	
	#Work out genetic lengths of unbroekn region
	rec.lengths.l <- gmap[position] - gmap[poss.bks.l] + 1e-8;
	
	#Work out amount of recombination in the interval with break taking edge cases into account
	if (poss.bks.l[1] == 1) {
		del.rec.l <- rep(0, l.b.l);
		del.rec.l[1] <- 1e-8;
		del.rec.l[2:l.b.l] <- gmap[poss.bks.l[2:l.b.l]] - gmap[poss.bks.l[2:l.b.l]-1] + 1e-8;
	} else {
		del.rec.l <- gmap[poss.bks.l] - gmap[poss.bks.l-1] + 1e-8;
	}
	
	llk.surf<-rep(0, length(t.grid));
	for (i in 1:length(t.grid)) {

		
		#tmp1 <-  -1 * t.grid[i] * rec.lengths.r + log(1 - exp(-1 * t.grid[i] * del.rec.r))  + p.bk.r +  (-1 * theta.per.site * t.grid[i] * lengths.r) + s.r * log(lengths.r * theta.per.site * t.grid[i]);
		#tmp2 <-  -1 * t.grid[i] * rec.lengths.l + log(1 - exp(-1 * t.grid[i] * del.rec.l))  + p.bk.l +  (-1 * theta.per.site * t.grid[i] * lengths.l) + s.l * log(lengths.l * theta.per.site * t.grid[i]);
		
		tmp1 <-  (s.r * log(t.grid[i] * theta.per.site * lengths.r) - (t.grid[i] * theta.per.site * lengths.r))   +   (log(1 - exp(-1 * t.grid[i] * del.rec.r)) - (t.grid[i] * rec.lengths.r))  +   p.bk.r
		tmp2 <-  (s.l * log(t.grid[i] * theta.per.site * lengths.l) - (t.grid[i] * theta.per.site * lengths.l))   +   (log(1 - exp(-1 * t.grid[i] * del.rec.l)) - (t.grid[i] * rec.lengths.l))  +   p.bk.l
		
		mx1<-max(tmp1);
		mx2<-max(tmp2);

		llk.surf[i] <- mx1 + mx2 + log(sum(exp(tmp1 - mx1))) + log(sum(exp(tmp2 - mx2)));
		
	}
	
	if (include.tij.prior) llk.surf<-llk.surf - t.grid;
	
	return(llk.surf);
	
}


#################
#End of Scripts
#################




