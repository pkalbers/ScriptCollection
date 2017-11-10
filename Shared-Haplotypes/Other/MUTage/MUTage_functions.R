

#Functions for mutation age


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


