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

		# first convert the sparse_array into a sparse_vector
		x <- x %>% as_sparse_vector()

		if (length(dims) == 1){

			x <- x %>% as_sparse_array()

		}else{

			subs <- NULL
			y <- x@subs
			for (i in 1:(length(dims) - 1)){
				n <- ((y - 1L) %% dims[i] + 1L)	# make the modulo from 1 to dims[i]. 
				y <- floor((y - 1L) / prod(dims[i])) + 1L
				subs <- cbind(subs, n)
			}
			subs <- cbind(subs, y)

			new(
				'sparse_array',
				subs = subs,
				vals = x@vals,
				dims = dims,
				dimnames = dimnames
			)
		}
	}
)
