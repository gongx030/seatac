windows_weight <- function(x, a = 25, b = 3/4){
	w <- (x / a)^b
	w[w > 1] <- 1
	w
} # windows_weight
