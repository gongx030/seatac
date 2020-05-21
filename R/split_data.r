split_data <- function(n, batch_size){
	starts <- seq(1, n, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > n] <- n
	cbind(start = starts, end = ends)
}
