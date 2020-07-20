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


