#' add_genomic_features
#'
#' Add genomic features
#'
#' @param x a Vplots object
#' @param y a GRanges object
#' @param x_name rowData(x) field for storing the features
#' @param y_name mcols(y) field for accessing the features
#' @param ... other parameters
#'
#' @return a Vplots object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_genomic_features',
	signature(
		x = 'Vplots',
		y = 'GRanges'
	), 
	function(
		x, 
		y,
		x_name = NULL,
		y_name = NULL,
		...
	){

		stopifnot(!is.null(y_name))
		stopifnot(!is.null(mcols(y)[[y_name]]))
		stopifnot(!is.null(x_name))

		if (is(mcols(y)[[y_name]], 'numeric')){

			rowData(x)[[x_name]] <- coverage(y, weight = y_name)[granges(x)] %>% as.matrix()

		}else if (is(mcols(y)[[y_name]], 'character')){

			y <- split(y, mcols(y)[[y_name]])
			y <- as(y, 'GRangesList')
			x <- add_genomic_features(x, y, x_name = x_name)
		}else
			stop(sprintf('unknown class type of mcols(y)[[%s]]', y_name))

		x

	}
)


#' add_genomic_features
#'
#' Add genomic features
#'
#' @param x a Vplots object
#' @param y a GRangesList object
#' @param x_name rowData(x) field for storing the features
#' @param ... other parameters
#'
#' @return a Vplots object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'add_genomic_features',
	signature(
		x = 'Vplots',
		y = 'GRangesList'
	), 
	function(
		x, 
		y,
		x_name = NULL,
		...
	){

		stopifnot(!is.null(x_name))

		n <- length(y)
		s <- names(y)
		h <- rep(1:n, sapply(y, length))
		y <- unlist(y)
		gr <- granges(x) %>%
			slidingWindows(width = x@bin_size, step = x@bin_size) %>%
			unlist()

		j <- matrix(1:(x@n_bins_per_window * length(x)), nrow = length(x), ncol = x@n_bins_per_window, byrow = TRUE) %>% c()
		gr <- gr[j]

		mm <- findOverlaps(gr, y) %>% as.matrix()
		mm <- sparseMatrix(i = mm[, 1], j = h[mm[, 2]], dims = c(length(gr), n), dimnames = list(NULL, s))
		dim(mm) <- c(length(x), x@n_bins_per_window * n)
		rowData(x)[[x_name]]<- mm
		x
	}
)
