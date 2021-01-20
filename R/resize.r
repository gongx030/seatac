#' Resize Vplots
#'
#' Resize a Vplots object
#'
#' @param x an input Vplots object
#' @param width Width of the sub-window
#' @param fix A character vector containing the values "start", "end", or "center" denoting what to use as an anchor for each element in x.
#' @return a new Vplots object
#' @export
#'
setMethod(
	'resize',
	signature(
		x = 'Vplots'
	),
	function(x, width, fix = 'center'){

		stopifnot(width %% x@bin_size == 0)

		if (width > x@window_size)
			return(x)

		if (fix == 'center'){
			center <- round(x@n_bins_per_window / 2)
			start <- round((x@window_size / x@bin_size - width / x@bin_size) / 2)
			end <- start + n_bins_per_block - 1
		}else
			stop(sprintf('unknown fix: %s', fix))

		n_bins_per_block <- as.integer(width / x@bin_size)

		d <- assays(x)$counts %>%
			summary() %>%
			as.matrix()

		hm <- matrix(1:x@n_bins_per_window, x@n_bins_per_window, x@n_intervals) 
		vm <- matrix(1:x@n_intervals, x@n_bins_per_window, x@n_intervals, byrow = TRUE)

		h <- hm[d[, 2]]
		v <- vm[d[, 2]]

		include <- h >= start & h <= end

		m <- matrix(1:(n_bins_per_block * x@n_intervals), n_bins_per_block, x@n_intervals)

		counts <- sparseMatrix(
			i = d[include, 1],
			j = m[cbind(h[include] - start + 1, v[include])],
			x = d[include, 3],
			dims = c(length(x), n_bins_per_block * x@n_intervals)
		)

		gr <- granges(x) %>%
			resize(fix = fix, width = width)

		se <- SummarizedExperiment(assays = list(counts = counts))
		SummarizedExperiment::rowRanges(se) <- gr
		rowData(se) <- rowData(x)

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
