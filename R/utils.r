#' add_track
#'
#' @param file a bigwig file
#'
add_track <- function(x, file, label){
	cvg <- rtracklayer::import(file, which = reduce(x), as = 'RleList')
  G <- sparseMatrix(
		i = 1:x@window_size,
		j = rep(1:x@n_bins_per_window, each = x@bin_size),
		x = rep(1 / x@bin_size, x@window_size),
		dims = c(x@window_size, x@n_bins_per_window)
	)
	y <- cvg[x] %>% as.matrix()
	y <- (y %*% G) %>% as.matrix()  # average signal in each genomic bin
	mcols(x)[[label]] <- y
	x
} # add_track

