setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))
setOldClass('kerastools.model.RModel')
setOldClass('tf_dataset')
setOldClass('tensorflow.tensor')
setOldClass('tensorflow.python.framework.sparse_tensor.SparseTensor')


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
		positions = 'numeric',
		n_samples = 'integer',
		samples = 'character'
	),
	contains = 'RangedSummarizedExperiment'
)


#' SummarizedVplots
#'
#' @export
setClass(
	'SummarizedVplots',
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
	contains = 'SummarizedExperiment'
)


setClass('Model', slot = c(model = 'kerastools.model.RModel'))
setClass('VaeModel', contains = 'Model')
setClass('SummarizedVaeModel', contains = 'Model')

#' VplotsList
#'
#' @export
#'
setClass(
	'VplotsList',
	contains = 'SimpleList'
)

#' SummarizedVplotsList
#'
#' @export
#'
setClass(
	'SummarizedVplotsList',
	contains = 'SimpleList'
)
