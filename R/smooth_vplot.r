
setMethod(
	'smooth_vplot',
	signature(
		x = 'GRanges'
	),
	function(x, theta = 1, ...){

		Z <- do.call('rbind', bplapply(1:length(x), function(i){
			z <- x[i]$counts %>%
				matrix(metadata(x)$n_bins_per_window, metadata(x)$n_intervals) %>%
				image.smooth(theta = theta) %>%
				pluck('z')

			z[z < 1e-7] <- 0

			z <- z / sum(z)
			c(z)
		}))

		Z <- as(Z, 'dgCMatrix')
			
		x$counts <- Z
		x

	}
)

