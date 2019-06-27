impute <- function(x, ...) UseMethod('impute')


impute.seatac_model <- function(x, filenames, which, num_windows_per_block = 1000, min_reads_per_window = 20){

	validate_bam(filenames)
  num_samples <- length(filenames)

	# dividing into small non-overlapping blocks; otherwise when making prediction keras will run out ouf resource on k40
	block_size <- x$window_size * num_windows_per_block
  blocks <- slidingWindows(which, width = block_size, step = block_size)
	blocks <- Reduce('c', blocks)

	bins <- NULL

	for (i in 1:length(blocks)){

		# the full bin sets in current block
		bs <- slidingWindows(blocks[i], width = x$bin_size, step = x$bin_size)
		bs <- unlist(bs)

		for (s in 1:num_samples){

			win <- getFragmentSizeMatrix(filenames[s], blocks[i], x$genome, x$window_size, x$bin_size, x$bin_size, x$fragment_size_range, x$fragment_size_interval, min_reads_per_window)

			if (!is.null(win)){

				sample_dim <- length(win)	# number of windows in current block
				flog.info(sprintf('prediction | block=%d/%d | region=%s:%d-%d | sample=%d/%d | sample_dim=%d', i, length(blocks), seqnames(blocks[i]), start(blocks[i]), end(blocks[i]), s, num_samples, sample_dim))

				# latent representation
				Z <- x$encoder(mcols(win)$counts %>% k_expand_dims() %>% tf$cast(tf$float32))$loc
				P <- x$latent_prior_model(NULL)$components_distribution$log_prob(Z %>% tf$reshape(shape(sample_dim, 1, x$latent_dim))) %>% as.matrix()

				# find the cluster membership of each window
				cls <- max.col(P)
				WC <- sparseMatrix(i = 1:length(win), j = cls, dims = c(length(win), x$n_components))	# window ~ cluster

				# predicted counts
				Xp <- x$decoder(Z)$distribution$probs %>% tf$reshape(shape(sample_dim, x$input_dim, x$feature_dim)) %>% as.array()
				Xp <- aperm(Xp, c(2, 1, 3))
				dim(Xp) <- c(sample_dim * metadata(win)$n_bins_per_window, metadata(win)$n_intervals)
 		 		bi <- slidingWindows(win, width = metadata(win)$bin_size, step = metadata(win)$bin_size)
				bi <- unlist(bi)
				g <- width(win) / metadata(win)$bin_siz	# # of bins per window

				BB <- as.matrix(findOverlaps(bs, bi)) # full bin set ~ covered bin set
				BB <- sparseMatrix(i = BB[, 1], j = BB[, 2], dims = c(length(bs), length(bi)))

				w <- 1 / rowSums(BB)	# number of coerved bins for each unique bin (since the windows may overlap)
				w[is.infinite(w)] <- 0
				BB <- Diagonal(x = w) %*% BB	# a weighted average mapping from covered bin to unique bin
				Xp <- as(BB %*% Xp, 'dgCMatrix')

				# observed counts per window
				Xi <- mcols(win)$counts
				Xi <- aperm(Xi, c(2, 1, 3))
				dim(Xi) <- c(sample_dim * metadata(win)$n_bins_per_window, metadata(win)$n_intervals)
				Xi <- as(BB %*% Xi, 'dgCMatrix')

				# the cluster membership
				# one bin can have multiple membership since the membership is assigned based on windows
				BW <- as.matrix(findOverlaps(bs, win)) # full bin set ~ window
				BW <- sparseMatrix(i = BW[, 1], j = BW[, 2], dims = c(length(bs), length(win)))
				BC <- BW %*% WC	# full bin set ~ cluster

				mcols(bs)$counts <- Xi
				mcols(bs)$predicted_counts <- Xp
				mcols(bs)$cluster <- BC
				bs <- bs[rowSums(BC) > 0]
				mcols(bs)$groups <- s

				if (is.null(bins))
					bins <- bs
				else
					bins <- c(bins, bs)
			}
		}
	}

	bins

} # predict.seatac_model

