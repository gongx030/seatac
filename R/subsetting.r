setMethod(
	'[', 
	'Vplots',
	function(x, i, j, ..., drop = FALSE){

		if (missing(i) && missing(j))
			return(x)

		if (!missing(j)){
			stop('subsetting by non-first dimension of a Vplots object is not implemented yet')
		}

		if (!missing(i)){
			x@dimdata[['grange']] <- x@dimdata[['grange']][i, , drop = FALSE]
			# need to update the id
		}

		callNextMethod(x, i, j, ..., drop = drop)

	}
)

setReplaceMethod(
	'[',
	'Vplots',
	function(x, i, j, ..., value){
		stop('not implemented yet')
	}
)
	
