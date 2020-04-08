setGeneric('get_motif_vplot', function(x, motifs, ...) standardGeneric('get_motif_vplot'))

#' get_motif_vplot
#'
#' @export
#'
setMethod(
	'get_motif_vplot',
	signature(
		x = 'GRanges',
		motifs = 'GRanges'
	),
	function(x, motifs, width = 320, ...){

		window_size <- width(x)

		if (length(unique(window_size)) > 1)
			stop('the window size of input data must be equal')

		window_size <- window_size[1]

		n_bins_per_motif <- width / metadata(x)$bin_size

		bins <- x %>%
			slidingWindows(width = metadata(x)$bin_size, step = metadata(x)$bin_size) %>%
			unlist()

		MB <- findOverlaps(resize(motifs, fix = 'center', width = 1), bins) %>%	# motif center ~ bins
			as.matrix()

		# there may be overlap among the peaks, so one motif center could be mapped to multiple peaks
		# if this happens, just randomly choose one motif/peak to proceed
		MB <- MB[!duplicated(MB[, 1]), ]

		MB <- MB[rep(1:nrow(MB), each = n_bins_per_motif), ]
		offsets <- 1:n_bins_per_motif - round(mean(n_bins_per_motif / 2))
		MB[, 2] <- MB[, 2] + offsets

		BF <- as.matrix(x$counts) 
		dim(BF) <- c(length(x), metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)	# batch size ~ bins per window ~ intervals ~ samples

		BF <- aperm(BF, c(2, 1, 3, 4))	# n_bins_per_window ~ batch size ~ n_intervals
		dim(BF) <- c(length(x) * metadata(x)$n_bins_per_window, metadata(x)$n_intervals, metadata(x)$n_samples)

		BF <- BF[MB[, 2], , ]	# length(motifs) * n_bins_per_motif ~ n_intervals ~ n_samples
		dim(BF) <- c(n_bins_per_motif, length(motifs), metadata(x)$n_intervals, metadata(x)$n_samples)
		BF <- aperm(BF, c(2, 1, 3, 4))	# length(motifs) ~ n_bins_per_motif ~ n_intervals ~ n_samples

		dim(BF) <- c(length(motifs), n_bins_per_motif * metadata(x)$n_intervals * metadata(x)$n_samples)
		motifs$counts <- as(BF, 'dgCMatrix')

		metadata(motifs)$n_samples <- metadata(x)$n_samples
		metadata(motifs)$samples <- metadata(x)$samples
		metadata(motifs)$fragment_size_range  <- metadata(x)$fragment_size_range
		metadata(motifs)$fragment_size_interval <- metadata(x)$fragment_size_interval
		metadata(motifs)$bin_size <- metadata(x)$bin_size
		metadata(motifs)$window_size <- width
		metadata(motifs)$n_intervals <- metadata(x)$n_intervals
		metadata(motifs)$n_bins_per_window <- n_bins_per_motif
		metadata(motifs)$breaks <- metadata(x)$breaks
		metadata(motifs)$centers <- metadata(x)$centers
		metadata(motifs)$n_motifs <- nlevels(motifs$tf_name)
		metadata(motifs)$motifs <- levels(motifs$tf_name)

		motifs

	}
) # get_motif_vplot


#' get_motif_vplot
#'
#' @export
#'
setMethod(
	'get_motif_vplot',
	signature(
		x = 'GRanges',
		motifs = 'GRangesList'
	),
	function(x, motifs, width = 320, ...){

		n <- sapply(motifs, length)
		motif_names <- names(motifs)
		motifs <- Reduce('c', motifs)
		motifs$tf_name <- factor(rep(motif_names, n), motif_names)
		motifs <- resize(motifs, width = width,  fix = 'center')

		get_motif_vplot(x, motifs, width, ...)
	}
)

