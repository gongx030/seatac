#' vplot
#'
#' @export
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

#' vplot
#'
#' @export
setMethod(
	'vplot',
	signature(
		x = 'SummarizedVplots'
	),
	function(
		x,
		field = 'counts',
		...
	){
		vplot_core(x, field, ...)
	}
)


vplot_core <- function(x, field, ...){
	z <- assays(x)[[field]]
	w <- 1 / rowSums(z)
	w[is.infinite(w)] <- 0
	z <- Diagonal(x = w) %*% z
	z <- colSums(z)
	z <- z / sum(z)
	z <- matrix(z, x@n_bins_per_window, x@n_intervals) 

	image(x = x@positions, y = x@centers, z, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), yaxt = 'n', xlab = '', ylab = 'fragment size', ...)
	abline(v = 0, lty = 2, col = 'yellow')
	axis(2, c(0, 100, 180, 247))
} # vplot_core
