setMethod(
	'as_vplot_autoencoder_rge_model',
	signature(
		x = 'vplot_autoencoder_model'
	),
	function(
		x,
		num_clusters,
		sigma,
		gamma,
		lambda
	){

		flog.info('converting vplot_autoencoder_model to vplot_autoencoder_rge_model')

		model <- x@data %>%
			build_model(
				name = 'vplot_autoencoder_rge_model',
				latent_dim = x@latent_dim,
				num_clusters = num_clusters,
				sigma = sigma,
				gamma = gamma,
				lambda = lambda
			)

		# initialize the weights
		k_random_uniform(c(1L, model@window_dim, model@interval_dim, 1L)) %>%
		model@encoder() %>%
		model@decoder()

		tfile <- tempfile(fileext = '.h5')

		x@encoder$save_weights(tfile)
		model@encoder$load_weights(tfile)

		x@decoder$save_weights(tfile)
		model@decoder$load_weights(tfile)

		model@data <- x@data
		model
	}
)

