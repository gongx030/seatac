setOldClass('tensorflow.tensor')
setOldClass('rpytools.call.VaeModel')
setOldClass('python.builtin.VaeModel')
setOldClass('tf_dataset')

setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

#' VaeModel
#'
#' @export
setClassUnion('VaeModel', members = c('rpytools.call.VaeModel', 'python.builtin.VaeModel'))


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
	contains = 'GRanges'
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


