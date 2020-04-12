setGeneric('as_sparse_array', function(x, ...) standardGeneric('as_sparse_array'))

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

