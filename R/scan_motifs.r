#' scan_motifs
#'
scan_motifs <- function(model, gr, chunk_size = 2^12, batch_size = 256, theta0 = 0.25, genome){

	window_dim <- length(gr)
  window_size <- metadata(gr)$window_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	starts <- seq(1, window_dim, by = chunk_size)
	ends <- starts + chunk_size - 1
	ends[ends > window_dim] <- window_dim
	n_chunk <- length(starts)

	bins <- NULL
	N <- window_size - model$kernel_size + 1

	for (b in 1:n_chunk){

		flog.info(sprintf('%d/%d', b, n_chunk))
		i <- starts[b]:ends[b]

		sequence <- unlist(strsplit(as.character(mcols(gr[i])$sequence), '')) %>%
			factor(c('A', 'C', 'G', 'T')) %>%
			as.numeric() %>%
			matrix(nrow = length(i), ncol = window_size, byrow = TRUE)
		
    sequence <- sequence - 1

		# a matrix indicating which fragment in which sequence passed the filter
		M <- model$sequence_motifs %>% predict(sequence, batch_size = batch_size, verbose = 1)
		M <- aperm(M, c(2, 1, 3))
		dim(M) <- c(prod(dim(M)[1:2]), dim(M)[3])
		j <- rowSums(M > theta0) > 0	# kernel that is activated in at least one filter
		M[M <= theta0] <- 0
		flog.info(sprintf('found %d active motifs', sum(j)))

		if (any(j)){
			bj <- gr[i] %>% slidingWindows(width = model$kernel_size, step = model$conv_strides) %>% unlist()
			bj <- bj[j]
			mcols(bj)$window_id <- rep(i, each = N)[j]
			mcols(bj)$bin_id <- rep(1:N, length(i))[j]
			mcols(bj)$motifs <- as(M[j, ], 'ngCMatrix')
			if (is.null(bins))
				bins <- bj
			else
				bins <- c(bins, bj)
		}
	}
	mcols(bins)$sequence <- getSeq(genome, bins)
	bins

} # scan_motifs

