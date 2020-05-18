#' vplot_predictor
#'
#' @export
#'
#' @author Wuming Gong
#'
setMethod(
	'vplot_predictor',
	signature(
		x = 'GRanges'
	),
	function(
		x,
		latent_dim = 25L,
		...
	){

		model <- new(
			'vplot_predictor_model',
			encoder = encoder_model(
				sample_embedding_dim = 20L,
				sequence_embedding_dim = 100L,
				latent_dim = latent_dim,
				sample_dim = as.integer(metadata(x)$n_samples)
			),
			decoder = decoder_model(
				window_dim = metadata(x)$n_bins_per_window,
				interval_dim = metadata(x)$n_intervals,
				filters0 = 64,
				filters = c(32L, 32L, 1L),
				kernel_size = c(3L, 3L, 3L),
				window_strides = c(2L, 2L, 2L),
				interval_strides = c(2L, 2L, 2L),
			),
			n_samples = as.integer(metadata(x)$n_samples),
			window_dim = as.integer(metadata(x)$n_bins_per_window),
			interval_dim = as.integer(metadata(x)$n_intervals),
			latent_dim = as.integer(latent_dim)
		)

		model %>% fit(x, ...)
		model

	}
)

