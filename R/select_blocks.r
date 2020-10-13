#' select_blocks
#' 
#' Prepare tf_dataset for model fitting
#' 
#' @return a `tf_dataset` object
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
select_blocks <- function(x, batch_size = 128L, min_reads = 0, n_blocks = NA, ...){

	batches <- cut_data(nrow(x), batch_size)

	if (!is.na(n_blocks))
		n_blocks_per_batch <- table(factor(sample(1:length(batches), n_blocks, replace = TRUE), 1:length(batches)))

	res <- list()
	for (j in seq_len(length(batches))){

		if (!is.na(n_blocks) && n_blocks_per_batch[j] == 0)
			next

		h <- batches[[j]]
		y <- x[h] %>% select_blocks_batch(...)

		# select blocks with minimum reads
		if (min_reads > 0 && !is.null(y$vplots)){
			include <- y$n >= min_reads
			y <- lapply(y, function(r) tf$boolean_mask(r, include))
		}

		# randomly sample non-empty blocks
		if (!is.na(n_blocks)){
			ns <- y[[1]]$shape[[1]] # n total blocks in current batch
			if (n_blocks_per_batch[j] < ns){
				include <- tf$cast(seq_len(ns) %in% sample(ns, n_blocks_per_batch[j]), tf$bool)
				y <- lapply(y, function(r) tf$boolean_mask(r, include))
			}
		}
		res[[j]] <- y
	}

	empty <- sapply(res, is.null)
	res <- res[!empty]
	fields <- names(res[[1]])
	res <- lapply(fields, function(f){
		tf$concat(lapply(res, function(xx) xx[[f]]), axis = 0L)
	})
	names(res) <- fields
	res

} # 


#' select_blocks_batch
#'
select_blocks_batch <- function(x, block_size, with_vplots = TRUE, with_kmers = FALSE, with_predicted_counts = FALSE, types = NULL){

	res <- list()

	if (with_kmers && is.null(x@kmers))
		stop('x$kmers must be available if with_kmers=TRUE')

	n_bins_per_block <- as.integer(block_size / x@bin_size)
	n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)
	n_blocks <- as.integer(n_blocks_per_window * nrow(x))

	if (with_vplots){

		y <- assays(x)$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32) %>%
				extract_blocks_from_vplot(n_bins_per_block) 
										
			y <- y %>%
				tf$reshape(c(y$shape[[1]] * y$shape[[2]], y$shape[[3]], y$shape[[4]], 1L))

			res$vplots <- y

			res$n <- y %>%
				tf$math$count_nonzero(c(1L, 3L), dtype = tf$float32) %>% # number of non-zero pixel per bar
				tf$math$count_nonzero(1L) %>% # number of bar that have non-zero pixels
				tf$cast(tf$float32)

			w <- y %>% tf$reduce_max(shape(1L), keepdims = TRUE)	# max reads per v-plot
			y <- y / tf$where(w > 0, w, tf$ones_like(w))	# scale so that the sum of each bar is one

			# add weight for each genomic bin
			w <- tf$reduce_sum(res$vplots, 1L, keepdims = TRUE) > 0
			w <- w %>% tf$cast(tf$float32)
			res$weight <- w
	}

	if (with_kmers){
		g <- rowData(x)$kmers %>%
			tf$cast(tf$int32) %>%
			tf$expand_dims(-1L) %>%
			tf$expand_dims(-1L) %>%
			tf$image$extract_patches(
				sizes = c(1L, block_size, 1L, 1L),
				strides = c(1L, x@bin_size, 1L, 1L),
				rates = c(1L, 1L, 1L, 1L),
				padding = 'VALID'
			) %>%
			tf$squeeze(axis = 2L)

		g <- g %>%
			tf$reshape(c(g$shape[[1]] * g$shape[[2]], block_size))

		res$kmers <- g
	}

	if (!is.null(types)){
		for (type in types){
			if (is.null(rowData(x)[[type]])){
				flog.warn(sprintf('field x$%s does not exist', type))
			}else{
				if (!is.matrix(rowData(x)[[type]])){
					flog.warn(sprintf('x$%s is not in a matrix format', type))
				}else{
					if (x@n_bins_per_window != ncol(rowData(x)[[type]])){
						flog.warn(sprintf('ncol(x$%s) must be equal to x@n_bins_per_window', type))
					}else{
						v <- rowData(x)[[type]] %>%
							tf$cast(tf$float32) %>%
							tf$expand_dims(-1L) %>%
							tf$expand_dims(-1L) %>%
							tf$image$extract_patches(
								sizes = c(1L, n_bins_per_block, 1L, 1L),
								strides = c(1L, 1L, 1L, 1L),
								rates = c(1L, 1L, 1L, 1L),
								padding = 'VALID'
							) %>%
							tf$squeeze(axis = 2L)

							v <- v %>%
								tf$reshape(c(v$shape[[1]] * v$shape[[2]], n_bins_per_block))
							v <- scale01(v)
							res[[type]] <- v
					}
				}
			}
		}
	}
	res
} 
