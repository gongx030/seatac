# add_kmers_counts
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_kmers_counts',
	signature(
		x = 'VplotsKmers',
		y = 'SummarizedKmers'
	),
	function(x, y, threshold = 1e-3, pseudo_count = 1){

		browser()

		WKP <- tf$SparseTensor(
			indices = cbind(
				rep(1:length(x), each = x@window_size) - 1, # batch (window)
				c(t(rowData(x)$kmers)),	# kmer
				rep(1:x@window_size, length(x)) - 1 # position (bin)
			),
			values = rep(1, length(x) * x@window_size),
			dense_shape = c(length(x), length(x@kmers), x@window_size)
		) # batch ~ kmer index ~ position

		Z <- assays(y)$counts %>%
			tf$cast(tf$float32)

		Z <- Z + pseudo_count
		N <-  tf$reduce_sum(Z, keepdims = TRUE)	# total counts
		p <- tf$reduce_sum(Z, axis = 0L, keepdims = TRUE) / N# probability that a voxel is positive
		Z <- (Z - tf$reduce_sum(Z, 1L, keepdims = TRUE) * p) / tf$math$sqrt(tf$reduce_sum(Z, 1L, keepdims = TRUE) * p * (1 - p))

		cutoff <- qnorm(1 - threshold)
		A <- tf$where(Z > 1, 1, 0)

											
											
#											c(t(rowData(x)$kmers)), 0:(length(x) * x@n_bins_per_window - 1)),
#			dense_shape = c(length(x@kmers), as.integer(length(x) * x@n_blocks_per_window))
#		)

		KB <- tf$SparseTensor(
		) # kmers ~ bins

		      KB <- tf$sparse$reorder(KB)
	}
)
