
#' show-Vplots
#'
setMethod(
	'show',
	signature(
		object = 'Vplots'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## n_intervals:%d\n', object@n_intervals))

		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))
		cat(sprintf('## n_bins_per_window:%d\n', object@n_bins_per_window))
	}
)


#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsFittedBlocksCluster'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## clusters: %s\n', object@k))
	}
)

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
