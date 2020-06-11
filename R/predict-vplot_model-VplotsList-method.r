#' vplot_parametric_vae_model
#'
setMethod(
	'predict',
	signature(
		model = 'vplot_model',
		x = 'VplotsList'
	),
	function(
		model,
		x,
		...
	){
		x_unlist <- unlist(x)
		x_unlist <- model %>% predict(x_unlist, ...)
		x <- relist(x_unlist, x)
		x@elementType <- 'Vplots'
		class(x) <- 'VplotsList'
		x

	}
)
