#' as_sparse_array
#'
#' Convert a dgCMatrix object into an sparse_array object
#'
#' @param x a dgCMatrix object
#'
setMethod(
	'as_sparse_array',
	signature(
		x = 'dgCMatrix'
	),
	function(x, ...){

		y <- summary(x) %>%
			as.matrix()
		new(
			'sparse_array',
			subs = y[, 1:2],
			vals = y[, 3],
			dims = dim(x),
			dimnames = dimnames(x)
		)
	}
)


#' as_sparse_array
#'
setMethod(
	'array_reshape',
	signature(
		x = 'sparse_vector'
	),
	function(x, ...){
		stop('not implemented yet')
	}
)

