#
# make look up tables for compressed haplotype data
#


number2binary = function(number, noBits) {
	binary_vector = rev(as.numeric(intToBits(number)))
	if(missing(noBits)) {
		return(binary_vector)
	} else {
		binary_vector[-(1:(length(binary_vector) - noBits))]
	}
}

pair2dos <- function(pair) {
	if (identical(pair, c(0,0))) return("-1")
	if (identical(pair, c(0,1))) return(" 0")
	if (identical(pair, c(1,0))) return(" 1")
	if (identical(pair, c(1,1))) return(" 2")
}

bits2four <- function(bits) {
	four <- split(bits, factor(c("a", "a", 
								 "b", "b", 
								 "c", "c", 
								 "d", "d")))
	
	dos <- c()
	for (pair in four) {
		dos <- c(dos, pair2dos(pair))
	}
	
	return(dos)
}


void <- function() number2binary(0, 2) # empty
dos0 <- function() number2binary(1, 2) # allele dosage = 0
dos1 <- function() number2binary(2, 2) # allele dosage = 1
dos2 <- function() number2binary(3, 2) # allele dosage = 2


for (i in 0:255) {
	bits <- number2binary(i, 8)
	four <- bits2four(bits)

	#cat(sprintf("{ %s }, ", paste(four, collapse=",")))
	
	#cat(sprintf("%s, ", (if(any(four == " 2")) " true" else "false" ) ))
	
	cat(sprintf("%d, ", length(which(four != "-1")) ))
	
	if ((i+1) %% 16 == 0) cat("\n")
}


