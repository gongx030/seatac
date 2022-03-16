#' Binding Vplots objects column-wise
#'
#' @param ... Vplots objects
#' @param deparse.level integer controlling the construction of labels in the case of non-matrix-like arguments
#' @return a Vplots object
#'
setMethod(
	'cbind',
	signature(
		'Vplots'
	), 
	function(
		...,
		deparse.level = 1	
	){

		objects <- list(...)

		ds <- do.call('rbind', lapply(objects, function(x) x@dimdata[['sample']]))
		ds[, 'id'] <- 1:nrow(ds)

		# bind the assays
		objects <- callNextMethod(...)
		objects@dimdata[['sample']] <- ds
		colData(objects)$sample <- rep(ds[, 'id'], each = nrow(objects@dimdata[['bin']])* nrow(objects@dimdata[['interval']]))
		objects
	}
)

