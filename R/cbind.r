#' cbind of multiple Vplots object
#'
#' @importFrom BiocGenerics cbind
#' @export
#'
setMethod(  
	'cbind',
	signature(
		'Vplots'
	),
  function(
		...,
		deparse.level = 1
	){
		browser()
	}
)

