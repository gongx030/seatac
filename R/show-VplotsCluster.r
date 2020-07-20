#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsCluster'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## clusters: %s\n', object@k))
	}
)
