
#' show-Vplots
#'
setMethod(
	'show',
	signature(
		object = 'Vplots'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## fragment_size_range:%d,%d\n', object@fragment_size_range[1], object@fragment_size_range[2]))
		cat(sprintf('## fragment_size_interval:%d\n', object@fragment_size_interval))
		cat(sprintf('## n_intervals:%d\n', object@n_intervals))

		cat(sprintf('## bin_size:%d\n', object@bin_size))
		cat(sprintf('## window_size:%d\n', object@window_size))
		cat(sprintf('## n_bins_per_window:%d\n', object@n_bins_per_window))
	}
)


#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsFitted'
	),
	function(object){
		callNextMethod()

		cat(sprintf('## model: %s\n', class(object@model)))
		cat(sprintf('## block_size: %s\n', object@model@block_size))
		cat(sprintf('## n_bins_per_block: %s\n', object@model@n_bins_per_block))
		cat(sprintf('## min_reads_per_block: %s\n', object@model@min_reads_per_block))
		cat(sprintf('## max_reads_per_pixel: %s\n', object@model@max_reads_per_pixel))
		cat(sprintf('## latent_dim: %s\n', object@model@latent_dim))
	}
)

#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsKmers'
	),
	function(object){
		callNextMethod()
		cat(sprintf('## kmers: %s\n', object@k))
	}
)

#' show-VplotsFitted
#'
setMethod(
	'show',
	signature(
		object = 'VplotsKmersFitted'
	),
	function(object){
		callNextMethod()
		cat(sprintf('## kmers: %s\n', object@k))
	}
)

