setMethod(
	'as_vplot_autoencoder_cluster_v2_model',
	signature(
		x = 'vplot_autoencoder_model'
	),
	function(
		x,
		sigma,
		gamma,
		num_clusters
	){

		flog.info('converting vplot_autoencoder_model to as_vplot_autoencoder_cluster_v2_model')

		slot_names <- slotNames(model)
		slot_names <- slot_names[!slot_names %in% c('encoder', 'decoder')]
		param <- lapply(slot_names, function(x) slot(model, x))
		names(param) <- slot_names

		param$sigma <- sigma
		param$gamma <- gamma
		param$num_clusters <- num_clusters

		model_new <- do.call(build_model, c(name = 'vplot_autoencoder_cluster_v2_model', x = NULL, param))

		# initialize the weights
		k_random_uniform(c(1L, model@n_bins_per_window, model@n_intervals, 1L)) %>%
			model_new@encoder() %>%
			model_new@decoder()

		tfile <- tempfile(fileext = '.h5')

		model@encoder$save_weights(tfile)
		model_new@encoder$load_weights(tfile)

		model@decoder$save_weights(tfile)
		model_new@decoder$load_weights(tfile)

		model_new
	}
)

