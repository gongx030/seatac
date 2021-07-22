#' vplot
#' 
#' Plot a Vplot object
#' 
#' @param x a Vplots object
#' @param field  The assays field that is used for Vplots (default: 'counts')
#' @param ... Arguments passed to vplot_core()
#'
#' @export
#'
setMethod(
	'vplot',
	signature(
		x = 'Vplots'
	),
	function(
		x,
		field = 'counts',
		...
	){
		vplot_core(x, field, ...)
	}
)


#' vplot_core
#' 
#' Plot a Vplot object
#' 
#' @param x a Vplots object
#' @param field  The assays field that is used for Vplots (default: 'counts')
#' @param ... Arguments passed to image()
#'
vplot_core <- function(x, field, ...){
	z <- assays(x)[[field]]
	w <- 1 / rowSums(z)
	w[is.infinite(w)] <- 0
	z <- Diagonal(x = w) %*% z
	z <- colSums(z)
	z <- z / sum(z)
	z <- matrix(z, x@n_bins_per_window, x@n_intervals) 

	image(x = x@positions, y = x@centers, z, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), xaxt = 'n', yaxt = 'n', xlab = '', ylab = 'fragment size', ...)
	abline(v = 0, lty = 2, col = 'yellow')
	axis(1, c(min(x@positions), 0, max(x@positions)), c(-x@window_size / 2, 0, x@window_size / 2))
	axis(2, seq(0, x@fragment_size_range[2], by = 100L))
} # vplot_core
