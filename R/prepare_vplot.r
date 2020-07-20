setMethod(
	'prepare_vplot',
	signature(
		x = 'Vplots'
	),
	function(
		x
	){

		y <- x$counts %>%
			as.matrix() %>%
			reticulate::array_reshape(c(    # convert into a C-style array
				length(x),
				x@n_intervals,
				x@n_bins_per_window,
				1L
			)) %>%
			tf$cast(tf$float32)

		y

	}
)
