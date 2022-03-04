#' dim-Vplots
#' 
#' Get the dimensions of a Vplots object
#
#' @param x a Vplots object
#'
#' @export
#'
setMethod(
	'dim',
	signature(
		x = 'Vplots'
	),
	function(x){
		sapply(x@dimdata, nrow)
	}
)


