setGeneric('read_fragment_size_profile', function(x, ga, ...) standardGeneric('read_fragment_size_profile'))

#' read_fragment_size_profile
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param ga  a GAlignments object.
#' @param max_isize fragment_size_range
#'
#' @export
#'

setMethod(
	'read_fragment_size_profile',
	signature(
		x = 'GRanges',
		ga = 'GAlignments'
	), 
	function(x, ga, max_isize = 1500){

		# compute the center point between PE reads
		# this is faster than using GAlignmentPairs
		ga <- ga[strand(ga) == '+']
		ga <- GRanges(
			seqnames = seqnames(ga), 
			range = IRanges(start(ga) + round(mcols(ga)$isize / 2), width = 1), 
			isize = mcols(ga)$isize
		)
		ga <- ga[ga$isize <= max_isize]

		WR <- as.matrix(findOverlaps(x, ga))	# window ~ reads
		WR <- sparseMatrix(i = WR[, 1], j = WR[, 2], dims = c(length(x), length(ga))) %>% as('dgCMatrix')	# window ~ reads
		RF <- sparseMatrix(i = 1:length(ga), j = ga$isize, dims = c(length(ga), max_isize)) %>% as('dgCMatrix')	# reads ~ fragment size
		mcols(x)$fragment_size_profile <- WR %*% RF # window ~ fragment size
		x
	}
) # read_fragment_size_profile

