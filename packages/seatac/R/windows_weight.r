windows_weight <- function(x, a = 25, b = 3/4){
	w <- (mcols(gr)$num_reads[b] / a)^b
	w[w > 1] <- 1
	w
} # windows_weight
