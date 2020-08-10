#' 
setMethod(
	'prepare_blocks',
	signature(
		x = 'Vplots',
		model = 'vplot_model'
	),
	function(
		 x,
		 model,
		 batch_size = 32L # v-plot per batch
	){

		starts <- seq(1, length(x), by = batch_size)
	  ends <- starts + batch_size - 1
		ends[ends > length(x)] <- length(x)
		n_batch <- length(starts)

		y <- lapply(seq_len(n_batch), function(i){

			b <- starts[i]:ends[i]

			xi <- x[b] %>%
				prepare_vplot() %>%
				extract_blocks_from_vplot(model@n_bins_per_block)

			xi <- xi %>%
				tf$reshape(c(xi$shape[0] * xi$shape[1], xi$shape[2], xi$shape[3], 1L))

			h <- xi %>%
				tf$math$count_nonzero(c(1L, 2L, 3L), dtype = tf$float32)

			xi %>% tf$boolean_mask(h > model@min_reads_per_block, axis = 0L)

		})

		y <- y %>% 
			tf$concat(0L) %>%
			tf$clip_by_value(clip_value_min = 0, clip_value_max = model@max_reads_per_pixel)

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
		 min_reads = NA
	){

		if (is.na(min_reads))
			min_reads <- model@min_reads_per_block

		y <- x %>%
			prepare_vplot() %>%
			extract_blocks_from_vplot(model@n_bins_per_block) 
		
		y <- y %>%
			tf$reshape(c(y$shape[0] * y$shape[1], y$shape[2], y$shape[3], 1L))

		h <- y %>%
			tf$math$count_nonzero(c(1L, 2L, 3L), dtype = tf$float32)

		include <- h >= min_reads

		y <- y %>% 
			tf$boolean_mask(include, axis = 0L) %>%
			tf$clip_by_value(clip_value_min = 0, clip_value_max = model@max_reads_per_pixel)

		g <- x$kmers %>%
			tf$cast(tf$int32) %>%
			tf$expand_dims(-1L) %>%
			tf$expand_dims(-1L) %>%
			tf$image$extract_patches(
				sizes = c(1L, model@n_bins_per_block * model@bin_size, 1L, 1L),
				strides = c(1L, model@bin_size, 1L, 1L),
				rates = c(1L, 1L, 1L, 1L),
				padding = 'VALID'
			) %>%
			tf$squeeze(axis = 2L)

		g <- g %>%
			tf$reshape(c(g$shape[0] * g$shape[1], model@block_size))

		g <- g %>% tf$boolean_mask(include, axis = 0L)

		list(vplots = y, kmers = g)

	}
)


