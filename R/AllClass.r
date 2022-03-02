setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))
setOldClass('kerastools.model.RModel')
setOldClass('tf_dataset')
setOldClass('tensorflow.tensor')
setOldClass('tensorflow.python.framework.sparse_tensor.SparseTensor')


#' Vplots
#'
#' @export
#'
setClass(
	'Vplots', 
	representation(
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
	contains = 'RangedSummarizedExperiment',
	prototype(
		n_samples = 1L
	)
)


#' Model
#'
setClass('Model', slot = c(model = 'kerastools.model.RModel'))


#' VplotsList
#'
#' @export
#'
setClass(
	'VplotsList',
	contains = 'SimpleList'
)

#' VaeModel
#'
setClass('VaeModel', contains = 'Model')

