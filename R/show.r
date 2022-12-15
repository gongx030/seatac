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

		sprintf('class: %s\n', class(object)) %>% cat()
		sprintf('dim: %s\n', paste(dim(object), collapse = ' ')) %>% cat()
		sprintf('assays(%d): %s\n', length(assays(object, withDimnames = FALSE)), paste0(names(assays(object, withDimnames = FALSE)), collapse = ' ')) %>% cat()
		sprintf('rowData(%d): %s\n', ncol(rowData(object)), paste0(colnames(rowData(object)), collapse = ' ')) %>% cat()
		sprintf('dimdata:\n') %>% cat()
		for (i in 1:length(object@dimdata)){
			sprintf('%s: names(%d): %s\n', names(object@dimdata)[i], ncol(object@dimdata[[i]]), paste(colnames(object@dimdata[[i]]), collapse = ' ')) %>% cat()
		}

		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))
	}
)

#' show-VplotsList
#' 
#' Print a VplotsList object
#
#' @param object a VplotsList object
#'
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'VplotsList'
	),
	function(object){

		cat('class: VplotsList\n')
		cat(sprintf('## Vplots objects: %d\n', length(object)))
		cat(sprintf('## samples: %d\n', object@n_samples))
		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))

	}
)
