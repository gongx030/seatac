setMethod(
	'as_vplot_autoencoder_cluster_model',
	signature(
		x = 'vplot_autoencoder_model'
	),
	function(
		x,
		num_clusters,
		sigma,
		gamma
	){

		flog.info('converting vplot_autoencoder_model to vplot_autoencoder_cluster_model')

		browser()

		build_model(
			name = 'vplot_autoencoder_cluster_model', 
			latent_dim = x@latent_dim,
			num_clusters = num_clusters,
			sigma = sigma,
			gamma = gamma
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
