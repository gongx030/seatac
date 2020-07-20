#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsFitted'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## model: %s\n', class(object@model)))
	}
)
