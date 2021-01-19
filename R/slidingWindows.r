#' slidingWindows
#'
#' Get the sliding V-plots
#'
#' @param x a Vplots object
#' @param width The width of the sliding window
#' @param step The step for each sliding window
#' @param batch_size The batch size of extraction operation (default: 4096L)
#' 
#' @return a new Vplot object which window_size is the specified width. The rowData and colData 
#' in the input Vplots object will be missing in the new Vplot object.
#'
#' @export
#'
setMethod(
	'slidingWindows',
  signature(
		x = 'Vplots'
	),
	function(x, width, step, batch_size = 4096L){

		stopifnot(width %% x@bin_size == 0)

		stopifnot(step %% x@bin_size == 0)

		if (width > x@window_size)
			stop(sprintf('width(%d) must be equal to or smaller than x@window_size(%d).', width, x@window_size))

		if (is.null(assays(x)$counts))
			stop('assays(x)$counts does not exist.')

		n_bins_per_block <- as.integer(width / x@bin_size)

		n_bins_per_step <- as.integer(step / x@bin_size)

		block_starts <- seq(1, x@n_bins_per_window - n_bins_per_block + 1, by = n_bins_per_step)
		block_ends <- block_starts + n_bins_per_block - 1
		n_blocks_per_window <- length(block_starts)

		batches <- cut_data(length(x), batch_size)
		counts <- list()

		for (i in 1:length(batches)){

			message(sprintf('slidingWindows | batch=%4.d/%4.d', i, length(batches)))

			b <- batches[[i]]

			counts[[i]] <- assays(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				tf$image$extract_patches(
					sizes = c(1L, x@n_intervals, n_bins_per_block, 1L),
					strides = c(1L, 1L, n_bins_per_step, 1L),
					rates = c(1L, 1L, 1L, 1L),
					padding = 'VALID'
				) %>%
				tf$squeeze(axis = 1L) %>%
				tf$reshape(c(length(b) * n_blocks_per_window, x@n_intervals, n_bins_per_block)) %>%
				tf$expand_dims(-1L) %>%
				tf$reshape(shape(length(b) * n_blocks_per_window, -1L)) %>%
				as.matrix() %>%
				as('dgCMatrix')
		}

		counts <- do.call('rbind', counts)

		gr <- granges(x) %>% 
			slidingWindows(width = width, step = step) %>% 
			unlist()

		se <- SummarizedExperiment(assays = list(counts = counts))
		SummarizedExperiment::rowRanges(se) <- gr
		rowData(se) <- rowData(x)[rep(1:length(x), each = n_blocks_per_window), ]

		new(
			'Vplots',
			se,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = block_size,
			n_intervals = as.integer(x@n_intervals),
			n_bins_per_window = n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, block_size, by = x@bin_size) - (block_size / 2),
			n_samples = x@n_samples,
			samples = x@samples
		)
	}
)

