#' read_fragment_size_profile
#'
#' @param x a GAlignments object.
#' @param windows a GRange object that define a set of genomic regions.
#' @param fragment_size_range fragment_size_range
#'
#' @export
#'
read_fragment_size_profile <- function(
	x, 
	windows, 
	max_isize = 1500
){

	breaks <- fragment_size_range[1]:fragment_size_range[2]

	# compute the center point between PE reads
	# this is faster than using GAlignmentPairs
	x <- x[strand(x) == '+']
	x <- GRanges(
		seqnames = seqnames(x), 
		range = IRanges(start(x) + round(mcols(x)$isize / 2), width = 1), 
		isize = mcols(x)$isize
	)
	x <- x[x$isize <= max_isize]

	WR <- as.matrix(findOverlaps(windows, x))	# window ~ reads
	WR <- sparseMatrix(i = WR[, 1], j = WR[, 2], dims = c(length(windows), length(x))) %>% as('dgCMatrix')	# window ~ reads
	RF <- sparseMatrix(i = 1:length(x), j = x$isize, dims = c(length(x), max_isize)) %>% as('dgCMatrix')	# reads ~ fragment size
	mcols(windows)$fragment_size_profile <- WR %*% RF # window ~ fragment size
	windows

} # read_fragment_size_profile

