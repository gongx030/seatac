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

		flog.info(sprintf('reading %s', model_file))
		model <- readRDS(model_file)

		if (class(model) == 'vplot_vae_model'){
			networks <- c('encoder', 'decoder', 'prior')
		}else if (class(model) == 'vplot_cvae_model'){
			networks <- c('encoder', 'decoder', 'prior', 'embedder')
		}else
			stop(sprintf('unknown model class: %s', class(model)))

		slot_names <- slotNames(model)
		slot_names <- slot_names[!slot_names %in% networks]
		param <- lapply(slot_names, function(x) slot(model, x))
		names(param) <- slot_names

		model <- do.call(build_model, c(name = class(model), x = NULL, param))


		# initialize the weights
		if (class(model) == 'vplot_vae_model'){

			z <- k_random_uniform(c(1L, model@n_intervals, model@n_bins_per_block, 1L)) %>%
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
