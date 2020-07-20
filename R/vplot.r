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

		z <- mcols(x)[[field]]
		z <- Diagonal(x = 1 / rowSums(z)) %*% z
		z <- colSums(z)
		z <- z / sum(z)
		z <- matrix(z, x@n_bins_per_window, x@n_intervals) 

		image(x = x@positions, y = x@centers, z, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), yaxt = 'n', xlab = '', ylab = 'fragment size', ...)
		abline(v = 0, lty = 2, col = 'yellow')
		axis(2, c(0, 100, 180, 247))
	}
)

