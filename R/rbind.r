setMethod(
	'rbind',
	signature(
		... = 'Vplots'
	),
	function(...){
		n_samples <- Reduce('+', sapply(list(...), function(x) x@n_samples))
		samples <- sapply(list(...), function(x) x@samples)
		x <- callNextMethod(...)
		x@n_samples <- n_samples
		x@samples <- samples
		x
	}
)
