setGeneric('array_reshape', function(x, dims, ...) standardGeneric('array_reshape'))

#' array_reshape
#'
setMethod(
	'array_reshape',
	signature(
		x = 'dgCMatrix',
		dims = 'numeric'
	),
	function(x, dims, dimnames = NULL, ...){

		if (prod(dims) != prod(dim(x)))
			stop(sprintf('dims [product %d] do not match the length of object [%d]', prod(dims),  prod(dim(x))))

		y <- x %>% 
			as_sparse_array() %>%
			array_reshape(dims)

	}
)

#' array_reshape
#'
setMethod(
	'array_reshape',
	signature(
		x = 'sparse_array',
		dims = 'numeric'
	),
	function(x, dims, dimnames = NULL, ...){

		if (prod(dims) != prod(dim(x)))
			stop(sprintf('dims [product %d] do not match the length of object [%d]', prod(dims),  prod(dim(x))))

		x <- x %>% as_sparse_vector()

		subs <- do.call('cbind', lapply(1:length(dims), function(i){
			if (i == 1)
				s <- rep(1:dims[i], times = prod(dims[-1]), each = 1)
			else if (i == length(dims))
				s <- rep(1:dims[i], times = 1, each = prod(dims[-length(dims)]))
			else
				s <- rep(1:dims[i], each = prod(dims[1:(i - 1)]), times = prod(dims[(i + 1):length(dims)]))
			s[x@subs]
		}))

		new(
			'sparse_array',
			subs = subs,
			vals = x@vals,
			dims = dims,
			dimnames = dimnames
		)
	}
)

