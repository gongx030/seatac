#' slidingWindows
#'
setMethod(
	'slidingWindows',
  signature(
		x = 'Vplots'
	),
	function(x, width, batch_size = 8192L){

		stopifnot(width %% x@bin_size == 0)

		width <- min(width, x@window_size)

		n_bins_per_block <- as.integer(width / x@bin_size)
		n_blocks_per_window <- as.integer(x@n_bins_per_window - n_bins_per_block + 1)
		n_blocks <- as.integer(n_blocks_per_window * length(x))

		if (!is.null(assays(x))){

			assays <- list()
			batches <- cut_data(length(x), batch_size)

			for (h in names(assays(x))){

				assays[[h]] <- list()

				for (i in 1:length(batches)){

					b <- batches[[i]]

					y <- assays(x[b])[[h]] %>%
						as.matrix() %>%
						reticulate::array_reshape(c(    # convert into a C-style array
							length(b),
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
				
					assays[[h]][[i]] <- y
				}

				assays[[h]] <- do.call('rbind', assays[[h]])
			}

			se <- SummarizedExperiment(assays = assays)

		}else{
			se <- SummarizedExperiment()
		}

		gr <- x@rowRanges %>% 
			slidingWindows(width = width, step = x@bin_size) %>% 
			unlist()

		SummarizedExperiment::rowRanges(se) <- gr

		# copy slots
		slots <- slotNames(x)
		slots <- slots[!slots %in% slotNames(se)]
		se <- new(class(x), se)
		for (s in slots){
			slot(se, s) <- slot(x, s)
		}
		se@window_size <- width
		se@n_bins_per_window <- n_bins_per_block
		se@positions <- seq(se@bin_size, se@window_size, by = se@bin_size) - (se@window_size / 2)

		rowdata_fields <- colnames(rowData(x))

		if (length(rowdata_fields) > 0){
			for (h in rowdata_fields){

				if (is.vector(rowData(x)[[h]])){

					y <- rep(rowData(x)[[h]], each = n_blocks_per_window)

				}else if (is.matrix(rowData(x)[[h]])){

					dim_h <- dim(rowData(x)[[h]])

					if (dim_h[2] == x@window_size){
						patch_width <- width
						patch_strides <- x@bin_size 
					}else if (dim_h[2] ==  x@n_bins_per_window){
						patch_width <- as.integer(width / x@bin_size)
						patch_strides <- 1L
					}else
						stop(sprintf('dim(rowData(x)$%s)[2] must be either %d or %d', h, x@window_size, x@n_bins_per_window))

					y <- rowData(x)[[h]] %>%
						as.matrix() %>%
						tf$cast(tf$float32) %>%
						tf$expand_dims(-1L) %>%
						tf$expand_dims(-1L) %>%
						tf$image$extract_patches(
							sizes = c(1L, patch_width, 1L, 1L),
							strides = c(1L, patch_strides, 1L, 1L),
							rates = c(1L, 1L, 1L, 1L),
							padding = 'VALID'
						) %>%
						tf$squeeze(axis = 2L)

					y <- y %>%
						tf$reshape(c(y$shape[[1]] * y$shape[[2]], patch_width)) %>%
						as.matrix()

					if (is(assays(x)[[h]], 'sparseMatrix'))
						y <- as(y, 'dgCMatrix')

				}else
					stop(sprintf('rowData(x)$%s must be a vector or a matrix', h))

				SummarizedExperiment::rowData(se)[[h]] <- y
			}
		}
		se
	}
)

