#' nol-Vplots
#' 
#' Get the number of columns of a Vplots object
#
#' @param x a Vplots object
#'
#' @export
#'
setMethod(
	'ncol',
	signature(
		x = 'Vplots'
	),
	function(x){
		prod(dim(x)[-1L])
	}
)


