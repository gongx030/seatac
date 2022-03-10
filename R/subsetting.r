#' Subsetting a Vplots object
#'
#' @param x a Vplots object
#' @param i row-wise indexing
#' @param j column-wise indexing
#' @param ... Additional arguments
#' @param drop Whether or not drop dimension
#'
setMethod(
	'[', 
	signature(
		x = 'Vplots', 
		i = 'ANY', 
		j = 'ANY'
	),
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

#' Subsetting a Vplots object
#'
#' @param x a Vplots object
#' @param i row-wise indexing
#' @param j column-wise indexing
#' @param ... Additional arguments
#' @param value The replacement value
#'
setReplaceMethod(
	'[',
	signature(
		x = 'Vplots', 
		i = 'ANY', 
		j = 'ANY',
		value = 'Vplots'
	),
	function(x, i, j, ..., value){
		stop('not implemented yet')
	}
)
	
