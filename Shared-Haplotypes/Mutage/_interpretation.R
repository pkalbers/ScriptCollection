

cat("\nReading data");
load("rea2.RData");
cat("\nDone");

#Get positions of mutations
POS<-real_dat[,1];

#Read in functions
#source("MUTage_functions.R")

#Set Ne
NE <- 10000  # N0


#Reorganise input files
HAP<-real_dat[,-1]; # d
rm(real_dat);

#Make genetic map
MAP<-approx(tav[,1],tav[,3], xout=POS, rule=2)$y*4*NE/100; # gmap

#Decide whether to include uncertainty about rec breaks
hard.r.bks<-F;


#Data params
marker.size <- nrow(HAP); # l
sample.size <- ncol(HAP); # n

select.k <- 100;  ##k<-40 # k

pairs.max <- 10000;

window.max <- 2000;  # win_max

#Grid for times (in 2Ne units)
times.seq <- seq(0, 5, length.out = 1001); # tt.seq
times.seq[1] <- times.seq[2] / 2;

times.pos <- rep(0, length(times.seq)); # tt.pos


AAC <- rowSums(HAP); # f


#Stuff for breakpoint function
tmp <- AAC / sample.size;

ff <- tmp^2 + (1-tmp)^2;  # cumulative map of log homozygosities
ff <- cumsum(log(ff));    

hh <- log(2*tmp*(1-tmp)); # per SNP list of log heterozygosities


#Mutation rate and N_e
mut.rate.per.gen<-1.2e-8;
n_e<-1e4;

theta.per.site <- 4 * n_e * mut.rate.per.gen;         #Ne is just used as an internal scaling factor

alpha_n <- sum(1 / (1:sample.size));

theta_est_Watterson <- marker.size / (alpha_n * diff(range(POS)));

theta_est_pwd <- sum(2 * AAC / sample.size * (1 - AAC / sample.size)) / diff(range(POS));

#Some storage space
OUT <- array(0, c(pairs.max, 10));  # op
colnames(OUT)<-c("ID1", "ID2", "Mut.1", "Mut.2", "S.L", "S.R", "BK.L", "BK.R", "last.S.L", "last.S.R");




###################




#Choose sites to analyse
focal.sites <- 50001:50001; # pp


EST <- array(0, c(length(focal.sites), 8));  # pt.estimate
colnames(EST) <- c("Post.mean", "Post.mode", "Post.median", "Post.2.5", "Post.97.5", "Robust", "Site", "Freq");

EST[,"Site"] <- POS[focal.sites];
EST[,"Freq"] <- AAC[focal.sites];




iterator <- 0;  # ct.var


#Now iterate over sites
for (focal.site in focal.sites) {
	
	aac = AAC[focal.site]
	
	if (aac > 0 && aac < sample.size) {
		
		iterator <- iterator + 1;
		
		cat("\nDoing site ", iterator, " of ", length(focal.sites), sep="");
		
		window.lower <- max(1, focal.site - window.max);  # win.low
		window.upper <- min(marker.size, focal.site + window.max);  # win.high
		
		#Work out who has mutation
		sharers <- which(HAP[focal.site, ] == 1);  # wi
		others <- setdiff(1:sample.size, sharers);  # wni
		### n.i <- length(sharers);
		
		#Sample pairs of haplotypes to analyse
		#Within clade
		# if (n.i > 1) {
		# 	tmp1 <- sample.pairs(sharers, select.k);
		# } else { 
		# 	tmp1 <- array(0, c(0,2));
		# }
		sharer.pairs <- sample.pairs(sharers, select.k);  # tmp1
		
		#Across clade
		other.pairs <- sample.pairs2(sharers, others, select.k);  # tmp2
		
		#Augment pairs with maximal sharing outside clade
		sharer.list <- sample.safe(sharer.pairs, select.k);  # who.to.analyse
		
		#l.i <- length(sharer.list);
		
		tmp3 <- array(0, c(length(sharer.list), 2));
		tmp3[,1] <- sharer.list;
		
		#   tmp3[,2]<-sample(wni, l.i, replace=TRUE);           #Currently just a random selection
		
		for (i in 1:length(sharer.list)) {
			
			if (focal.site < marker.size) {
				dp<-apply(HAP[(focal.site + 1) : window.upper, others, drop=F] != HAP[(focal.site + 1) : window.upper, tmp3[i,1]], 2, find.first.true.element);
			} else {
				dp<-rep(0, length(others));
			}
			
			if (focal.site > 1) {
				dm <- apply(HAP[(focal.site - 1) : window.lower, others, drop=F] != HAP[(focal.site - 1) : window.lower, tmp3[i,1]], 2, find.first.true.element);
			} else {
				dm<-rep(0, length(others));
			}
			
			tmp3[i,2] <- others[which.max(dp+dm)];
		}
		
		#Count pairs
		#n.k <- nrow(sharer.pairs) + nrow(other.pairs) + nrow(tmp3);
		#ii <- 1:n.k;
		
		len = nrow(sharer.pairs) + nrow(other.pairs) + nrow(tmp3)  # n.k
		rng = 1 : len  # ii
		
		#Start filling op
		OUT[rng, "ID1"] <- c(sharer.pairs[,1], other.pairs[,1], tmp3[,1]);
		OUT[rng, "ID2"] <- c(sharer.pairs[,2], other.pairs[,2], tmp3[,2]);
		
		#Randomise
		#for (i in 1:n.k) {
		#op[i,1:2]<-sample(n, 2, replace=F);
		#}
		
		OUT[rng, "Mut.1"] <- as.vector(HAP[focal.site, OUT[rng, "ID1"]]);
		OUT[rng, "Mut.2"] <- as.vector(HAP[focal.site, OUT[rng, "ID2"]]);
	
		###OUT[rng, 5:10] <- 0;
		
		
		
		
		#Now go through each pair getting breaks and number of mutations
		for (i in 1:len) {
			
			tmp <- which(HAP[, OUT[i, "ID1"]] != HAP[, OUT[i, "ID2"]]);
			
			###f.focal<-f[p];  ### using aac
			
			
			#Analyse to RHS of pos
			breakpoints <- tmp[tmp > focal.site];  # b.r
			
			flag <- TRUE;
			
			if (length(breakpoints) == 0) {
				flag <- FALSE; 
				break.r <- marker.size + 1;
			}
			
			j <- 0;
			count <- 0;  # ct
			last.site <- focal.site;  # last.s
			
			while(flag) {
				j <- j+1;
				
				aac.sharers <- sum(HAP[breakpoints[j], sharers]);  # f.within
				
				if ((aac.sharers < aac) && (aac.sharers < AAC[breakpoints[j]])) {
					flag<-FALSE;
					break.r <- breakpoints[j];
				} else {
					count <- count+1; 
					last.site <- breakpoints[j];
				}
				
				if (j==length(breakpoints) && flag) {
					flag <- FALSE; 
					break.r <- marker.size + 1;
				}
			}
			
			OUT[i, "BK.R"] <- break.r;	#This is the location of the first difference - i.e. the point at which the break is detected
			
			OUT[i, "S.R"] <- count;
			
			OUT[i, "last.S.R"] <- last.site;
			
			if (OUT[i, "Mut.2"] == 0) {
				OUT[i, "S.R"] <- OUT[i, "S.R"]+1;	#Add +1 for focal SNP if discordant
			}
			
			
			
			#Analyse to LHS of pos
			breakpoints <- tmp[tmp < focal.site];  # b.l
			
			flag<-TRUE;
			
			if (length(breakpoints)==0) {
				flag<-FALSE; 
				break.l <- 0;
			}
			
			j <- length(breakpoints) + 1;
			count <- 0;
			last.site <- focal.site;
			
			while(flag) {
				j<-j-1;
				
				aac.sharers <- sum(HAP[breakpoints[j], sharers]);
				
				if ((aac.sharers < aac) && (aac.sharers < AAC[breakpoints[j]])) {
					flag<-FALSE;
					break.l <- breakpoints[j];
				} else {
					count <- count + 1; 
					last.site <- breakpoints[j];
				}
				
				if (j==1 && flag) {
					flag <- FALSE; 
					break.l <- 0;
				}
			}
			
			OUT[i, "BK.L"] <- break.l;	#This is the location of the first difference - i.e. the point at which the break is detected
			
			OUT[i, "S.L"] <- count;
			
			OUT[i, "last.S.L"] <- last.site;
		}
		
		
		
		
		#
		#
		#Now infer age distribution
		#
		#
		
		if (hard.r.bks) {
			A_PAR <- 1 + OUT[rng, "S.L"] + OUT[rng, "S.R"] + as.integer(OUT[rng, "BK.L"] > 0) + as.integer(OUT[rng, "BK.R"] <= marker.size);  # a.params
			#    a.params<-1+op[ii,5]+op[ii,6];
			
			bkr <- OUT[rng, "BK.R"];  # tmp1
			bkr[bkr > marker.size] <- marker.size;
			
			bkl <- OUT[rng, "BK.L"];  # tmp2
			bkl[bkl < 1] <- 1;
			
			mut.length <- (POS[bkr] - POS[bkl] + 1) * theta.per.site;
			rec.length <- (MAP[bkr] - MAP[bkl]);
			
			B_PAR <- 1 + mut.length + rec.length;  # b.params
			#    b.params<-1+mut.length;
		}
		
		
		times.pos[]<-0;
		
		plot(times.seq, rep(0, length(times.seq)), type="n", ylim=c(0,1), xlab="Time", ylab="CCF", log="x", main=paste("Pos = ", POS[focal.site], ": f = ", AAC[focal.site], sep=""));
	
		v.low<-0;
		v.high<-0;
		
		
		for (i in 1:len) {
			
			OUT_I <- OUT[i,];  # op.i
			position <- focal.site;
			
			if (!hard.r.bks) {
				llk.surf.i <- llk.tmrca.pair(OUT[i,], 
																	 position = focal.site, 
																	 t.grid = times.seq, 
																	 gmap = MAP, 
																	 pos = POS, 
																	 ff = ff, 
																	 hh = hh, 
																	 theta.per.site = theta.per.site, 
																	 debug.flag=FALSE);
				
				mx<-max(llk.surf.i);
				tmp1<-exp(llk.surf.i-mx);
				
				tmp1<-tmp1/sum(tmp1);
				cdf.i<-cumsum(tmp1);
				cdf.i<-1e-8+(1-2e-8)*cdf.i
			}
			
			if (OUT[i,4]==1) {
				if (hard.r.bks) {
					times.pos<-times.pos+pgamma(times.seq, A_PAR[i], B_PAR[i], log=T, lower=T);
					lines(times.seq, pgamma(times.seq, A_PAR[i], B_PAR[i], lower=T), col="blue");
					
					pg = pgamma(times.seq, A_PAR[i], B_PAR[i], lower=T)
					v.low<-v.low+log(times.seq[which.min(abs(pg-0.5))]);
					
				} else {
					times.pos<-times.pos+log(cdf.i);
					
					lines(times.seq, cdf.i, col="blue");
					
					v.low<-v.low+log(times.seq[which.min(abs(cdf.i-0.5))]);
				}
			}
			if (OUT[i,4]==0) {
				if (hard.r.bks) {
					times.pos<-times.pos+pgamma(times.seq, A_PAR[i], B_PAR[i], log=T, lower=F);
					lines(times.seq, pgamma(times.seq, A_PAR[i], B_PAR[i], lower=F), col="red");
					
					pg = pgamma(times.seq, A_PAR[i], B_PAR[i], lower=F)
					v.high<-v.high+log(times.seq[which.min(abs(pg-0.5))]);
					
				} else {
					times.pos<-times.pos+log(1-cdf.i);
					
					lines(times.seq, 1-cdf.i, col="red");
					
					v.high<-v.high+log(times.seq[which.min(abs(cdf.i-0.5))]);
				}
			}
			if (sum(is.na(times.pos))>0) break();
		}
		
		mx<-max(times.pos);
		times.pos<-exp(times.pos-mx);
		times.pos<-times.pos/sum(times.pos);
		
		#plot(times.seq, times.pos, log="x", xlab="Age", ylab="Relative Posterior", type="l", 
		#    col="blue", main=paste("Pos = ", POS[focal.site], ": f = ", AAC[focal.site], sep=""));
		lines(times.seq, times.pos/max(times.pos), type="l", col="black", lwd=2);
		#abline(v=true.age[p], lty="dotted", col="red");
		
		
		n.within<-sum(OUT[rng, "Mut.2"]);
		
		if (n.within>0) {
			v.low <- v.low / n.within;
		}
		
		v.high <- v.high / (max(len - n.within, 1));
		
		if (n.within > 0) {
			con.est <- (v.low + v.high) / 2;
		} else {
			con.est <- v.high - log(2);
		}
		
		abline(v=exp(c(v.low, v.high, con.est)), col=c("blue", "red", "lightblue"));
		
		EST[iterator, 6]<-con.est;
		
		
		#Get pseudo-95% ETPI
		times.cum<-cumsum(times.pos);  # tt.cum
		oo<-approx(x=c(0,times.cum), y=c(0,times.seq), xout=c(0.025, 0.975));
		
		EST[iterator, 1] <- sum(times.seq * times.pos);
		EST[iterator, 2] <- times.seq[which.max(times.pos)];
		EST[iterator, 3] <- times.seq[which.min(abs(cumsum(times.pos)-0.5))];
		
		EST[iterator, 4:5] <- oo$y;  
		
	}
	
}



write.table(EST, file="output.txt", quote=F, row=F, col=T);



