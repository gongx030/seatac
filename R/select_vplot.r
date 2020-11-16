#'
#' @export
#'
setMethod(
	'select_vplot',
	signature(
		x = 'Vplots'
	),
	function(x, assays = NULL, fields = NULL){

		assay_names <- names(assays(x))
		assay_names <- assay_names[assay_names %in% assays]
		SummarizedExperiment::assays(x) <- assays(x)[assay_names]

		field_names <- colnames(rowData(x))
		field_names <- field_names[field_names %in% fields]
		SummarizedExperiment::rowData(x) <- rowData(x)[, field_names, drop = FALSE]

		x
	}
)
