setMethod(
	'load_model',
	signature(
		dir = 'character'
	),
	function(
		dir,
		...
	){


		if (!file.exists(dir))
			stop(sprintf('%s does not exist dir'))

		model_file <- sprintf('%s/model.rds', dir)

		if (!file.exists(model_file))
			stop(sprintf('%s does not exist', model_file))

		# initialize the weights
		if (class(model) == 'vplot_vae_model'){

			z <- tf$random$uniform(c(1L, model@n_intervals, model@n_bins_per_block, 1L)) %>%
				model@encoder() 
		
			z$sample() %>%
				model@decoder()

		}else if (class(model) == 'vplot_cvae_model'){

			x <- k_random_uniform(c(1L, model@n_intervals, model@n_bins_per_block, 1L)) 

			g <- k_zeros(c(1L, model@block_size), dtype = tf$int32)

			h <- g %>% 
				model@embedder()

			z <- list(x, h) %>%
				model@encoder()

			list(z$sample(), h) %>%
				model@decoder()
		}

		for (s in networks){
			s_file <- sprintf('%s/%s.h5', dir, s)
			slot(model, s)$load_weights(s_file)
			flog.info(sprintf('reading %s', s_file))
		}

		model	

	}
)
