#' c
#' 
#' Concatenate two or more Vplots objects
#' @param x a Vplots object
#'
#' @export
#'
#' @return a Vplots object
#'
setMethod(
	'c',
  signature(
		x = 'Vplots'
	),
	function(x, ...){

		assay_names <- names(assays(x))
		rowdata_names <- colnames(rowData(x))

		x <- list(x, ...)

		a <- lapply(assay_names, function(h) Reduce('rbind', lapply(x, function(xx) assays(xx)[[h]])))
		names(a) <- assay_names

		gr <- Reduce('c', lapply(x, function(xx) rowRanges(xx)))

		se <- SummarizedExperiment(assays = a)
		SummarizedExperiment::rowRanges(se) <- gr

		# copy slots
		slots <- slotNames(x[[1]])
		slots <- slots[!slots %in% slotNames(se)]
		se <- new(class(x[[1]]), se)
		for (s in slots){
			slot(se, s) <- slot(x[[1]], s)
		}
		se
	}
)
