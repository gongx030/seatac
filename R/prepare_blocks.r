#' 
setMethod(
	'prepare_blocks',
	signature(
		x = 'Vplots',
		model = 'vplot_vae_model'
	),
	function(
		 x,
		 model,
		 min_reads = 0
	){

		block_size <- model@encoder$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)

		y <- x %>%
			prepare_vplot() %>%
			extract_blocks_from_vplot(n_bins_per_block) 
		
		y <- y %>%
			tf$reshape(c(y$shape[0] * y$shape[1], y$shape[2], y$shape[3], 1L))

		h <- y %>%
			tf$math$count_nonzero(c(1L, 2L, 3L), dtype = tf$float32)

		include <- h >= min_reads

		y <- y %>% 
			tf$boolean_mask(include, axis = 0L) %>%
			tf$clip_by_value(clip_value_min = 0, clip_value_max = 100)

		y

	}
)


#' 
setMethod(
	'prepare_blocks',
	signature(
		x = 'VplotsKmers',
		model = 'vplot_cvae_model'
	),
	function(
		 x,
		 model,
		 min_reads = 0 
	){

		block_size <- model@encoder$block_size
		n_bins_per_block <- as.integer(block_size / x@bin_size)

		y <- x %>%
			prepare_vplot() %>%
			extract_blocks_from_vplot(n_bins_per_block) 
		
		y <- y %>%
			tf$reshape(c(y$shape[0] * y$shape[1], y$shape[2], y$shape[3], 1L))

		h <- y %>%
			tf$math$count_nonzero(c(1L, 2L, 3L), dtype = tf$float32)

		include <- h >= min_reads

		y <- y %>% 
			tf$boolean_mask(include, axis = 0L) %>%
			tf$clip_by_value(clip_value_min = 0, clip_value_max = 100)

		g <- x$kmers %>%
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
			tf$reshape(c(g$shape[0] * g$shape[1], block_size))

		g <- g %>% tf$boolean_mask(include, axis = 0L)

		list(vplots = y, kmers = g)

	}
)


