setGeneric('read_vplot', function(x, peaks, ...) standardGeneric('read_vplot'))

#' read_vplot
#'
#' @param x a GAlignments object.
#' @param peaks a GRange object that define a set of genomic regions.
#' @param window_size The size of each genomic window for training the model.
#' @param bin_size The bin size
#' @param fragment_size_range fragment_size_range
#' @param fragment_size_interval fragment_size_interval
#'
#' @export
#'
setMethod(
	'read_vplot',
	signature(
		x = 'GAlignments',
		peaks = 'GRanges'
	), 
	function(
		x, 
		peaks, 
		max_isize = 650,
		...
	){

		if (peaks %>% width() %>% unique() %>% length() > 1)
			stop('each peak must have the same width')

		w <- width(peaks)[1]	# assuming that all peaks have the same width

		# compute the center point between PE reads
		# this is faster than using GAlignmentPairs
		x <- x[strand(x) == '+']
		x <- GRanges(
			seqnames = seqnames(x), 
			range = IRanges(start(x) + round(mcols(x)$isize / 2), width = 1), 
			isize = mcols(x)$isize,
		)
		x <- x[x %over% peaks]
		x <- x[x$isize <= max_isize]

		positions <- peaks %>%
			slidingWindows(width = 1, step = 1) %>%
			unlist()	# 

		PR <- findOverlaps(positions, x) %>% as.matrix()
		PR <- sparseMatrix(i = PR[, 1], j = PR[, 2], dims = c(length(positions), length(x))) %>% as('dgCMatrix')	# positions ~ reads
		RF <- sparseMatrix(i = 1:length(x), j = mcols(x)$isize, dims = c(length(x), max_isize))	%>% as('dgCMatrix') # reads ~ fragment size
		FP <- t(PR %*% RF)	# fragment_size ~ positions
		dim(FP) <- c(w * max_isize, length(peaks))
		mcols(peaks)$counts <- t(FP)
		peaks
	}
)

