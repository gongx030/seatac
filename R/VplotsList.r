#' VplotsList
VplotsList <- function(...){

	listData <- list(...)

	if (length(listData) == 0L) {
		unlistData <- new('Vplots', GRanges())
	}else {
		if (!all(sapply(listData, is, 'Vplots')))
			stop("all elements in '...' must be Vplots objects")

		unlistData <- suppressWarnings(do.call('c', unname(listData)))
	}

	x <- relist(unlistData, listData)
	x@elementType <- 'Vplots'
	class(x) <- 'VplotsList'
	x

}
