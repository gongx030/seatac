#' slidingWindows
#'
setMethod(
	'slidingWindows',
  signature(
		x = 'Vplots'
	),
	function(x){

		block_size <- x@block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)
		n_blocks <- as.integer(n_blocks_per_window * nrow(x))

		if (!is.null(assays(x))){

			assays <- list()

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

		}else{
			se <- SummarizedExperiment()
		}

		gr <- x@rowRanges %>% 
			slidingWindows(width = block_size, step = x@bin_size) %>% 
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

		if (!is.null(colnames(rowData(x)))){

			for (h in colnames(rowData(x))){

				dim_h <- dim(rowData(x)[[h]])

				if (dim_h[2] == x@window_size){
					y <- rowData(x)[[h]] %>%
						as.matrix() %>%
						tf$cast(tf$float32) %>%
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

				}else if (dim_h[2] == n_blocks_per_window){

					if (length(dim(rowData(x)[[h]])) == 2)
						target_shape <- shape(dim_h[1] * dim_h[2], 1L)
					else
						target_shape <- shape(dim_h[1] * dim_h[2], 1L, dim_h[-c(1:2)])

					y <- rowData(x)[[h]] %>%
						tf$cast(tf$float32) %>%
						tf$reshape(target_shape) %>%
						as.array()

				}else
					stop(sprintf('dim(rowData(x)$%s)[2] must be either %d or %d', h, x@window_size, n_blocks_per_window))

				SummarizedExperiment::rowData(se)[[h]] <- y
			}
		}
		se
	}
)
