#' Summarize V-plot surround each kmers
#' 
#' Summarize the aggregated V-plot surround each kmers
#'
#' @param x a VplotsKmers object 
#' @param block_size the window size for the V-plot surrounding each k-mer
#' @param batch_size batch size 
#' @return a SummarizedKmers object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'summarize_kmers',
	signature(
		x = 'VplotsKmers'
	), 
	function(
		x, 
		block_size,
		batch_size = 512L
	){

		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)

		KV <- tf$zeros(shape(length(x@kmers), x@n_intervals, n_bins_per_block, 1L))	# kmer ~ Vplot
		n <- tf$zeros(shape(length(x@kmers)))

		batches <- cut_data(length(x), batch_size)
		valid <- (n_bins_per_block / 2):(x@n_bins_per_window - n_bins_per_block/ 2)	# blocks that are completely contained within the window

		for (i in 1:length(batches)){

			message(sprintf('summarize_kmers | batch=%6.d/%6.d', i, length(batches)))

			b <- batches[[i]]

			BV <- assays(x[b])$counts %>%
				as.matrix() %>%
				reticulate::array_reshape(c(    # convert into a C-style array
					length(b),
					x@n_intervals,
					x@n_bins_per_window,
					1L
				)) %>%
				tf$cast(tf$float32) %>%
				extract_blocks_from_vplot(n_bins_per_block)	%>% # batch ~ n_blocks_per_window ~ n_intervals ~ n_bins_per_block ~ 1L
				tf$reshape(c(length(b) * n_blocks_per_window, -1L))

			KB <- tf$SparseTensor(
				indices = cbind(c(t(rowData(x[b])$kmers[, valid])), 0:(length(b) * n_blocks_per_window - 1)),
				values = rep(1, length(b) * n_blocks_per_window), 
				dense_shape = c(length(x@kmers), as.integer(length(b) * n_blocks_per_window))
			) # kmers ~ bins
			KB <- tf$sparse$reorder(KB)

			KV_i <- tf$sparse$sparse_dense_matmul(KB, BV) %>%
				tf$reshape(c(length(x@kmers), x@n_intervals, n_bins_per_block, 1L))

			KV <- KV + KV_i
			n <- n + KB %>% tf$sparse$reduce_sum(1L)

		}

		counts <- KV %>%
			tf$reshape(c(length(x@kmers), -1L)) %>%
			as.matrix()
		rownames(counts) <- x@kmers

		se <- SummarizedExperiment(
			assays = list(counts = counts)
		)
		SummarizedExperiment::rowData(se)$n <- as.numeric(n)

		se <- new(
			'SummarizedKmers',
			se,
			fragment_size_range  = x@fragment_size_range,
			fragment_size_interval = x@fragment_size_interval,
			bin_size = x@bin_size,
			window_size = block_size,
			n_intervals = x@n_intervals,
			n_bins_per_window = n_bins_per_block,
			breaks = x@breaks,
			centers = x@centers,
			positions = seq(x@bin_size, block_size, by = x@bin_size) - (block_size / 2),
			kmers = x@kmers,
			k = x@k
		)

		se <- compute_z_score(se)
		se
	}
) # summarize_kmers


