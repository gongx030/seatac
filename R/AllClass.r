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
		dimdata = 'list'
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
	contains = 'SimpleList',
	representation(
		n_samples = 'integer',
		fragment_size_range  = 'integer',
		fragment_size_interval = 'integer',
		bin_size = 'integer',
		window_size = 'integer'
	),
	prototype(
		n_samples = 1L
	),
	validity = function(object){

		msg <- NULL

		window_sizes <- sapply(object, function(obj) obj@window_size)
		if (length(unique(window_sizes)) > 1){
			msg <- 'All Vplots objects should have identical window sizes'
		}

		bin_sizes <- sapply(object, function(obj) obj@bin_size)
		if (length(unique(bin_sizes)) > 1){
			msg <- 'All Vplots objects should have identical bin sizes'
		}

		fragment_size_intervals <- sapply(object, function(obj) obj@fragment_size_interval)
		if (length(unique(fragment_size_intervals)) > 1){
			msg <- 'All Vplots objects should have identical fragment size intervals'
		}

		fragment_size_starts <- sapply(object, function(obj) obj@fragment_size_range[1])
		if (length(unique(fragment_size_starts)) > 1){
			msg <- 'All Vplots objects should have identical fragment size starts'
		}

		fragment_size_ends <- sapply(object, function(obj) obj@fragment_size_range[2])
		if (length(unique(fragment_size_ends)) > 1){
			msg <- 'All Vplots objects should have identical fragment size ends'
		}

		return(msg)

	}
)

VplotsList <- function(...){

	objects <- new('VplotsList', ...)
	objects@n_samples <- sapply(objects, function(obj) nrow(obj@dimdata$sample)) %>% sum()
	objects@fragment_size_range <- objects[[1]]@fragment_size_range
	objects@fragment_size_interval <- objects[[1]]@fragment_size_interval
	objects@bin_size <- objects[[1]]@bin_size
	objects@window_size <- objects[[1]]@window_size
	objects

}



#' VaeModel
#'
setClass('VaeModel', contains = 'Model')

