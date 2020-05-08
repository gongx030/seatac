#' vplot_autoencoder
#'
#' @export
#'
#' @author Wuming Gong
#'
setMethod(
	'vplot_autoencoder',
	signature(
		x = 'GRanges'
	),
	function(
		x,
		latent_dim = 10L,
		epochs = 50,
		batch_size = 128,
		steps_per_epoch = 50,
		...
	){

		model <- new(
			'vplot_autoencoder_model',
			encoder = encoder_model(
				latent_dim = latent_dim,
				filters = c(32L, 32L, 32L),
				kernel_size = c(3L, 3L, 3L),
				window_strides = c(2L, 2L, 2L),
				interval_strides = c(2L, 2L, 1L)
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

		model %>% fit(x, batch_size = batch_size, epochs = epochs, steps_per_epoch = steps_per_epoch)

		z <- model %>% encode(x, batch_size = batch_size)

		metadata(x)$model <- model
		x$latent <- z

		x
	}
)

