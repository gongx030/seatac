setMethod(
	'add_nucleosome_signal',
	signature(
		x = 'Vplots'
	),
	function(x){

		if (is.null(x$predicted_counts))
			stop('predicted_counts field does not exist')

		nfr <- which(x@centers > 0 & x@centers <= 100)
		mono <- which(x@centers >= 180 & x@centers <= 247)

		y <- x$predicted_counts
		dim(y) <- c(length(x), x@n_bins_per_window, x@n_intervals)

		x$predicted_nfr <- y[, , nfr, drop = FALSE] %>% rowMeans(dims = 2)
		x$predicted_mono_nucleosome <- y[, , mono, drop = FALSE] %>% rowMeans(dims = 2)

		x$nucleosome_signal <- (x$predicted_mono_nucleosome + 1e-3) / (x$predicted_nfr + 1e-3)

		x$nucleosome_signal <- (x$nucleosome_signal - rowMins(x$nucleosome_signal)) / (rowMaxs(x$nucleosome_signal) - rowMins(x$nucleosome_signal))
		x

	}
)
