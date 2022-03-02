#' cbind of multiple Vplots object
#'
#' @importFrom BiocGenerics cbind
#' @importFrom MatrixGenerics rowRanges
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

		args <- list(...)

		x <- callNextMethod(...)
		browser()

	}
)

