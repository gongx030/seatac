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

#		callNextMethod()

		sprintf('class: %s', class(object)) %>% message()
		sprintf('dim: %d %d %d %d', nrow(object), x@n_samples, x@n_intervals, x@n_bins_per_window) %>% message()
		sprintf('assays(%d): %s', length(assays(object, withDimnames = FALSE)), paste0(names(assays(object, withDimnames = FALSE)), collapse = ' ')) %>% message()
		sprintf('dimdata:') %>% message()
		for (i in 1:length(object@dimdata)){
			sprintf('%s: names(%d): %s', names(object@dimdata)[i], ncol(object@dimdata[[i]]), paste(colnames(object@dimdata[[i]]), collapse = ' ')) %>% message()
		}

		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## n_intervals:%d\n', object@n_intervals))
		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))
		cat(sprintf('## n_bins_per_window:%d\n', object@n_bins_per_window))
	}
)


