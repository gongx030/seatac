#' Smooth the V-plot
#'
#' @param x GRanges object
#'
setMethod(
	'smooth_vplot',
	signature(
		x = 'GRanges'
	),
	function(x, theta = 1, ...){

		w <- 1 / rowSums(x$counts)
		w[is.infinite(w)] <- 0
		x$smoothed_counts <- Diagonal(x = w) %*% x$counts
		n <- rowSums(x$counts)

		Z <- do.call('rbind', bplapply(1:length(x), function(i){
			
			if (n[i] > 0){
				z <- x[i]$smoothed_counts %>%
					matrix(metadata(x)$n_bins_per_window, metadata(x)$n_intervals) %>%
					image.smooth(theta = theta) %>%
					pluck('z')

				z[z < 1e-7] <- 0

				z <- z / sum(z)

				c(z) %>%
					matrix(1, metadata(x)$n_bins_per_window * metadata(x)$n_intervals) %>%
					as('dgCMatrix')

			}else{

				x[i]$smoothed_counts

			}
		}))

		x$smoothed_counts <- Z
		x

	}
)

