setGeneric('compute_variability', function(x, ...) standardGeneric('compute_variability'))


setMethod(
	'compute_variability',
	signature(
		x = 'SummarizedExperiment'
	),
	function(x, ...){

		Z <- assays(x)$z
		Z[is.na(Z)] <- 0
		Z[is.infinite(Z)] <- 0
		dim(Z) <- c(metadata(x)$n_motifs, metadata(x)$n_bins_per_window,  metadata(x)$n_intervals, metadata(x)$n_samples)

		Z_mean <- rowMeans(Z, dims = 3)
		Z_mean <- array(
			c(Z_mean),
			dim = c(
				metadata(x)$n_motifs,
				metadata(x)$n_bins_per_window,
				metadata(x)$n_intervals,
				metadata(x)$n_samples
			)
		)

		Z_var <- rowSums((Z - Z_mean)^2, dims = 3) / metadata(x)$n_samples

		P <- pchisq(
			(metadata(x)$n_samples  - 1) * Z_var, 
			df = (metadata(x)$n_samples - 1), 
			lower.tail = FALSE
		)
										                  

		browser()

		-log10(P[n, , ]) %>% image(col = colorpanel(100, low = 'black', high = 'yellow'), breaks = c(seq(0, -log10(1e-3), length.out = 100), 10000))
																			                 
	}
)


