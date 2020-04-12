setGeneric('as_sparse_vector', function(x, ...) standardGeneric('as_sparse_vector'))

setMethod(
	'as_sparse_vector',
	signature(
		x = 'dgCMatrix'
	),
	function(x, ...){

		dim(x) <- c(prod(dim(x)), 1)
		x <- summary(x) %>%
			as.matrix()
		new(
			'sparse_vector',
			subs = x[, 1],
			vals = x[, 3]
		)
	}
)

setMethod(
	'as_sparse_vector',
	signature(
		x = 'sparse_array'
	),
	function(x, ...){

		p <- sapply(2:length(dim(x)), function(i) prod(dim(x)[1:(i - 1)]))
		subs <- rowSums(t(t(x@subs[, -1, drop = FALSE] - 1) * p)) + x@subs[, 1]
		new(
			'sparse_vector',
			subs = subs,
			vals = x@vals
		)
	}
)

