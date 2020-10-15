#' add_track
#'
#' @param file a bigwig file
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_track',
	signature(
		x = 'Vplots',
		object = 'character'
	),
	function(x, object, label){
		cvg <- rtracklayer::import(object, which = reduce(x@rowRanges), as = 'RleList')
	  G <- sparseMatrix(
			i = 1:x@window_size,
			j = rep(1:x@n_bins_per_window, each = x@bin_size),
			x = rep(1 / x@bin_size, x@window_size),
			dims = c(x@window_size, x@n_bins_per_window)
		)
		y <- cvg[x@rowRanges] %>% as.matrix()
		y <- (y %*% G) %>% as.matrix()  # average signal in each genomic bin
		SummarizedExperiment::rowData(x)[[label]] <- y
		x
	}
) # add_track

