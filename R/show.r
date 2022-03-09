#' show-Vplots
#' 
#' Print a Vplots object
#
#' @param object a Vplots object
#'
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'Vplots'
	),
	function(object){

		sprintf('class: %s', class(object)) %>% message()
		sprintf('dim: %s', paste(dim(object), collapse = ' ')) %>% message()
		sprintf('assays(%d): %s', length(assays(object, withDimnames = FALSE)), paste0(names(assays(object, withDimnames = FALSE)), collapse = ' ')) %>% message()
		sprintf('rowData(%d): %s', ncol(rowData(object)), paste0(colnames(rowData(object)), collapse = ' ')) %>% message()
		sprintf('dimdata:') %>% message()
		for (i in 1:length(object@dimdata)){
			sprintf('%s: names(%d): %s', names(object@dimdata)[i], ncol(object@dimdata[[i]]), paste(colnames(object@dimdata[[i]]), collapse = ' ')) %>% message()
		}

		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))
	}
)
