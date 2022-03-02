#' vplot
#' 
#' Plot a Vplot object
#' 
#' @param x a Vplots object
#' @param field  The assays field that is used for Vplots (default: 'counts')
#' @param ncol number of columns (default: 2L)
#'
#' @export
#'
setMethod(
	'vplot',
	signature(
		x = 'Vplots'
	),
	function(
		x,
		field = 'counts',
		ncol = 2L
	){
		vplot_core(x, field, ncol)
	}
)


#' vplot_core
#' 
#' Plot a Vplot object
#' 
#' @param x a Vplots object
#' @param field  The assays field that is used for Vplots (default: 'counts')
#' @param ncol number of columns (default: 2L)
#' @import ggplot2
#'
vplot_core <- function(x, field, ncol = 2){

	stopifnot(!is.null(assays(x)[[field]]))

	z <- assays(x)[[field]]

	df <- do.call('rbind', lapply(1:x@n_samples, function(i){
		zi <- z[, colData(x)$batch == x@samples[i], drop = FALSE]
		w <- 1 / rowSums(zi)
		w[is.infinite(w)] <- 0
		zi <- Diagonal(x = w) %*% zi
		zi <- colSums(zi)
		zi <- zi / sum(zi)
		data.frame(value = zi, bin = x@positions[colData(x)$bin], interval = x@centers[colData(x)$interval], batch = x@samples[i])
	}))	%>%
		mutate(batch = factor(batch, x@samples))

	df %>%
		ggplot(aes(x = bin, y = interval, fill = value)) +
			geom_raster() +
			facet_wrap( ~ batch, ncol = ncol) +
			theme(panel.spacing = unit(2, 'lines')) +
			geom_vline(xintercept = 0, linetype = 'dotted', color = 'yellow') +
			scale_fill_gradientn(colors = colorpanel(100, low = 'blue', mid = 'white', high = 'red')) +
			scale_x_continuous(
				breaks = c(-x@window_size / 2, 0, x@window_size / 2),
				limits = c(-x@window_size / 2, x@window_size / 2),
				expand = c(0, 0)
			) + 
			scale_y_continuous(
				breaks = seq(0,  x@fragment_size_range[2], by = 100L),
				limits = c(0, max(x@centers)),  
				expand = c(0, 0)
			) + 
			xlab('') + ylab('fragment size')

} # vplot_core
