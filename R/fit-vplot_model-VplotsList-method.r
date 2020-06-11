#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_model',
		x = 'VplotsList'
	),
	function(
		model,
		x,
		...
	){
		x <- unlist(x)
		fit(model, x, ...)
	}
)


