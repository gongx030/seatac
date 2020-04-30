#' get_motif_vplot
#'
#' @export
#'
setMethod(
	'get_motif_vplot',
	signature(
		x = 'GRanges',
		pwms = 'PWMatrixList'
	),
	function(x, pwms, width = 320, genome, min_reads = 20, max_reads = 500, ...){

		window_size <- width(x)

		if (length(unique(window_size)) > 1)
			stop('the window size of input data must be equal')

		window_size <- window_size[1]

		motifs <- matchMotifs(
			pwms,
			resize(x, fix = 'center', width = peak_size - width),
			genome = genome,
			out = 'positions'
		)

		n <- sapply(motifs, length)
		motif_names <- names(motifs)
		motifs <- Reduce('c', motifs)
		motifs$tf_name <- factor(rep(motif_names, n), motif_names)
		motifs <- resize(motifs, width = width,  fix = 'center')

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

		BF <- x$counts %>% as_sparse_array()
		BF <- BF %>% array_reshape(
			c(
				length(x),
				metadata(x)$n_bins_per_window,
				metadata(x)$n_intervals,
				metadata(x)$n_samples
			)
		)
		BF <- BF %>% 
			array_permute(c(2, 1, 3, 4))# n_bins_per_window ~ batch size ~ n_intervals

		BF <- BF %>% array_reshape(
			c(
				length(x) * metadata(x)$n_bins_per_window, 
				metadata(x)$n_intervals * metadata(x)$n_samples
			)
		)

		BF <- sparseMatrix(i = BF@subs[, 1], j = BF@subs[, 2], x = BF@vals, dims = BF@dims)

		BF <- BF[MB[, 2], ] # length(motifs) * n_bins_per_motif ~ n_intervals *  n_samples

		BF <- BF %>% as_sparse_array()

		BF <- BF %>% array_reshape(
			c(
				n_bins_per_motif, 
				length(motifs), 
				metadata(x)$n_intervals, 
				metadata(x)$n_samples
			)
		)

		BF <- BF %>% 
			array_permute(c(2, 1, 3, 4)) # length(motifs) ~ n_bins_per_motif ~ n_intervals ~ n_samples

		BF <- BF %>% array_reshape(
			c(
				length(motifs), 
				n_bins_per_motif * metadata(x)$n_intervals * metadata(x)$n_samples
			)
		) # length(motifs) ~ n_bins_per_motif * n_intervals * n_samples


		BF <- sparseMatrix(i = BF@subs[, 1], j = BF@subs[, 2], x = BF@vals, dims = BF@dims)

		motifs$counts <- BF

		motifs$n_reads <- rowSums(motifs$counts)
		motifs <- motifs[motifs$n_reads >= min_reads & motifs$n_reads < max_reads] # some region have high number of reads

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
		metadata(motifs)$positions <- seq(metadata(motifs)$bin_size, metadata(motifs)$window_size, by = metadata(motifs)$bin_size) - (metadata(motifs)$window_size / 2)

		motifs

	}
) # get_motif_vplot


