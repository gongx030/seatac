setGeneric('as_dgCMatrix', function(x, ...) standardGeneric('as_dgCMatrix'))

setMethod(
	'as_dgCMatrix',
	signature(
		x = 'sparse_array'
	),
	function(x){

		if (length(dim(x)) != 2)
			stop('dim must be 2')

		sparseMatrix(
			i = x@subs[, 1],
			j = x@subs[, 2],
			x = x@vals,
			dims = x@dims,
			dimnames = x@dimnames
		)
	}
)

