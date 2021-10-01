#' get_nucleosome_diff
#' 
#' Get nucleosome difference between two V-plot objects
#' 
#' @param x a Vplots object
#' @param y a Vplots object
#' @param fragment_size_threshold Fragment size threshold where the any fragment size above is considered as nucleosome reads and below is considered as nucleosome free reads
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
get_nucleosome_diff <- function(x, y, fragment_size_threshold = 150){

	stopifnot(!is.null(assays(x)$predicted_counts))
	stopifnot(!is.null(assays(y)$predicted_counts))

	is_nucleosome <- x@centers >= fragment_size_threshold

	v_from <- assays(x)$predicted_counts %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
		scale01()

	v_to <- assays(y)$predicted_counts %>%
		tf$cast(tf$float32) %>%
		tf$reshape(shape(-1L, x@n_intervals, x@n_bins_per_window, 1L)) %>%
		scale01()

	nuc_to <- v_to %>%
		tf$boolean_mask(is_nucleosome, 1L) %>%
		tf$reduce_sum(1L)

	nfr_to <- v_to %>%
		tf$boolean_mask(!is_nucleosome, 1L) %>%
		tf$reduce_sum(1L)

	nuc_from <- v_from %>%
		tf$boolean_mask(is_nucleosome, 1L) %>%
		tf$reduce_sum(1L)

	nfr_from <- v_from %>%
		tf$boolean_mask(!is_nucleosome, 1L) %>%
		tf$reduce_sum(1L)

	nuc_nfr_from <- log((nuc_from + 1e-3) / (nfr_from + 1e-3))
	nuc_nfr_to <- log((nuc_to + 1e-3) / (nfr_to + 1e-3))

	nuc_diff <- 2 / (1 + exp(-1 * (nuc_nfr_to - nuc_nfr_from))) - 1
	nuc_diff <- nuc_diff %>% tf$squeeze(2L)
	nuc_diff

}
