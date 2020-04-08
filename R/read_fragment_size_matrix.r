setGeneric('get_vplot', function(x, ga, ...) standardGeneric('get_vplot'))

#' get_vplot
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param ga a GAlignments object.
#' @param bin_size The bin size
#' @param fragment_size_range fragment_size_range
#' @param fragment_size_interval fragment_size_interval
#'
#' @export
#'
setMethod(
	'get_vplot',
	signature(
		x = 'GRanges',
		ga = 'GAlignmentsList'
	), 
	function(
		x, 
		ga, 
		bin_size = 5,
		fragment_size_range = c(50, 690), 
		fragment_size_interval = 10
	){

		window_size <- width(x)

		if (length(unique(window_size)) > 1)
			stop('the window size of input data must be equal')
		
		window_size <- window_size[1]

 		n_bins_per_window <- window_size / bin_size

		breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

		n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval

		wb <- cbind(rep(1:length(x), n_bins_per_window), 1:(length(x)* n_bins_per_window))	# windows ~ bins, sorted by window

		bins <- x %>%
			slidingWindows(width = bin_size, step = bin_size) %>%
			unlist()

		# compute the center point between PE reads
		# this is faster than using GAlignmentPairs

		mcols(x)$counts <- do.call('cbind', lapply(ga, function(g){	
		
			g <- g[strand(g) == '+']
			g <- GRanges(
				seqnames = seqnames(g), 
				range = IRanges(start(g) + round(mcols(g)$isize / 2), width = 1), 
				isize = mcols(g)$isize
			)

			g$fragment_size <- as.numeric(cut(g$isize, breaks))	# discretize the fragment size
			g <- g[!is.na(g$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"

			CF <- sparseMatrix(i = 1:length(g), j = g$fragment_size, dims = c(length(g), n_intervals))  # read center ~ fragment size
			BC <- as.matrix(findOverlaps(bins, g))	# bins ~ read center
			BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(g))), 'dgCMatrix') # bins ~ read center
			BF <- BC %*% CF  # bins ~ fragment size
			BF <- as.matrix(BF[wb[, 2], ])
			dim(BF) <- c(n_bins_per_window, length(x), n_intervals)	# convert BF into an array with n_bins_per_window ~ batch_size ~ n_intervals
			BF <- aperm(BF, c(2, 1, 3)) # batch_size, n_bins_per_window ~ n_intervals
			dim(BF) <- c(length(x), n_bins_per_window * n_intervals)
			as(BF, 'dgCMatrix')

		}))	# batch_size ~ n_bins_per_window * n_intervals * n_samples

		metadata(x)$n_samples <- length(ga)
		metadata(x)$samples <- names(ga)
		metadata(x)$fragment_size_range  <- fragment_size_range
		metadata(x)$fragment_size_interval <- fragment_size_interval
		metadata(x)$bin_size <- bin_size
		metadata(x)$window_size <- window_size 
		metadata(x)$n_intervals <- n_intervals
		metadata(x)$n_bins_per_window <- n_bins_per_window 
		metadata(x)$breaks <- breaks
		metadata(x)$centers <- centers

		x
	}

) # get_vplot

