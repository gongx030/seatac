#' slidingWindows
#'
setMethod(
	'slidingWindows',
  signature(
		x = 'Vplots'
	),
	function(x, width, step = 1L){

		block_size <- as.integer(width)
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)
		n_blocks <- as.integer(n_blocks_per_window * nrow(x))

		assays <- list()

		# for assays
		for (h in names(assays(x))){

			y <- assays(x)[[h]] %>%
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

			y <- y %>% 
				tf$reshape(c(y$shape[[1]], -1L)) %>% 
				as.matrix()

			if (is(assays(x)[[h]], 'sparseMatrix'))
				y <- as(y, 'dgCMatrix')

			assays[[h]] <- y
		}

		se <- SummarizedExperiment(assays = assays)

		gr <- x@rowRanges %>% 
			slidingWindows(width = block_size, step = step) %>% 
			unlist()

		SummarizedExperiment::rowRanges(se) <- gr

		# copy slots
		slots <- slotNames(x)
		slots <- slots[!slots %in% slotNames(se)]
		se <- new(class(x), se)
		for (s in slots){
			slot(se, s) <- slot(x, s)
		}
		se@window_size <- block_size
		se@n_bins_per_window <- n_bins_per_block
		se@positions <-  seq(se@bin_size, se@window_size, by = se@bin_size) - (se@window_size / 2)

		fields <- colnames(rowData(x))
		if (length(fields) > 0){
			for (h in fields){
				y <- rowData(x)[[h]] %>%
					as.matrix() %>%
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

				y <- y %>%
					tf$reshape(c(y$shape[[1]] * y$shape[[2]], block_size)) %>%
					as.matrix()

				if (is(assays(x)[[h]], 'sparseMatrix'))
					y <- as(y, 'dgCMatrix')

				SummarizedExperiment::rowData(se)[[h]] <- y
			}
		}
		se
	}
)

