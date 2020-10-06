setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))
setOldClass('kerastools.model.RModel')


#' Vplots
#'
#' @export
setClass(
	'Vplots', 
	slot = c(
		fragment_size_range  = 'integer',
		fragment_size_interval = 'integer',
		bin_size = 'integer',
		window_size = 'integer',
		n_intervals = 'integer',
		n_bins_per_window = 'integer',
		breaks = 'numeric',
		centers = 'numeric',
		positions = 'numeric'
	),
	contains = 'RangedSummarizedExperiment'
)

#' VplotsKmers
#'
#' @export
setClass(
	'VplotsKmers',
	slot = c(
		kmers = 'character',
		k = 'integer'
	),
	contains = 'Vplots'
)

setClass(
	'VaeModel',
	slot = c(
		model = 'kerastools.model.RModel'
	)
)

		
