#' array_permute 
#'
setMethod(
	'array_permute',
	signature(
		x = 'sparse_array',
		perm = 'numeric'
	),
	function(x, perm, ...){

		if (length(perm) != length(dim(x)))
			stop(sprintf('perm is of wrong length %d (!= %d)', length(perm), length(dim(x))))

		if (any(perm > length(dim(x))))
				stop('value out of range')

		new(
			'sparse_array',
			subs = x@subs[, perm],
			vals = x@vals,
			dims = x@dims[perm],
			dimnames = x@dimnames[perm]
		)
	}
)
