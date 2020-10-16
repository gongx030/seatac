#' sample_blocks 
#' 
#' Randomly sample genomic blocks
#'
#' @param x a Vplots object.
#' @param block_size size of the block.
#' @param batch_size number of windows in one batch. 
#' @param min_reads minimum read number per block
#' @param max_reads minimum read number per block
#' @param n_blocks number of blocks to return
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
sample_blocks <- function(x, n_blocks = 1, block_size, batch_size = 128L, min_reads = 0, max_reads = 100000){

	batches <- cut_data(nrow(x), batch_size)

	if (!is.na(n_blocks))
		n_blocks_per_batch <- table(factor(sample(1:length(batches), n_blocks, replace = TRUE), 1:length(batches)))

	res <- list()
	for (j in seq_len(length(batches))){

		h <- batches[[j]]
		y <- x[h] %>% 
			slidingWindows(width = block_size, step = x@bin_size)

		message(sprintf('sample_blocks | batch=%4.d/%4.d | n total=%7.d | n sampled=%7.d', j, length(batches), length(y), n_blocks_per_batch[j]))

		n <- rowSums(assays(y)$counts)	# number of reads per block
		y <- y[n >= min_reads & n <= max_reads]

		if (length(y) > 0){
			y <- sample(y, n_blocks_per_batch[j])
			res[[j]] <- y
		}
	}
	empty <- sapply(res, is.null)
	res <- res[!empty]

	Reduce('c', res)

} # 
