#' extract_blocks
#'
setMethod(
	'extract_blocks',
	signature(
		x = 'VplotsFitted'
	),
	function(x){

		g <- matrix(1:(x@n_bins_per_window * x@n_intervals), x@n_bins_per_window, x@n_intervals)
		u <- do.call('cbind', lapply(1:x@model@n_blocks_per_window, function(i) c(g[i:(i + x@model@n_bins_per_block - 1), ])) )
		u <- c(u)

		gr <- x %>%
			slidingWindows(width = x@model@block_size, step = x@bin_size) %>%
			unlist()

		y <- x$counts[, u]
		y <- t(y)
		dim(y) <- c(x@model@n_bins_per_block * x@n_intervals, length(x) * x@model@n_blocks_per_window)
		y <- t(y)
		gr$counts <- y

		y <- x$predicted_counts[, u]
		y <- t(y)
		dim(y) <- c(x@model@n_bins_per_block * x@n_intervals, length(x) * x@model@n_blocks_per_window)
		y <- t(y)
		gr$predicted_counts <- y

		z <- x$latent
    z <- aperm(z, c(2, 1, 3))
    dim(z) <- c(length(x) * x@model@n_blocks_per_window, x@model@latent_dim)

		gr$latent <- z

		new(
			'VplotsFittedBlocks',
			gr,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = x@model@block_size,
			n_intervals = x@n_intervals,
			n_bins_per_window = x@model@n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, x@model@block_size, by = x@bin_size) - (x@model@block_size / 2),
			model = x@model
		)
	}
)

