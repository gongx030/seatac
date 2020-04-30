#' get_counts
#' 
#' Get the read counts for genomic intervals
#'
#' @param x a GRange object that define a set of genomic regions.
#' 
#' @export
#'
setMethod(
	'get_counts',
	signature(
		x = 'GRanges',
		filenames = 'character',
		genome = 'BSgenome'
	),
	function(x, filenames, ...){

		if (is.null(names(filenames)))
			stop('filenames must be a named vector')

		if (any(duplicated(names(filenames))))
			stop('names(filenames) must be unique')

		browser()
		X <- do.call('cbind', bplapply(filenames, function(file){	
			g <- read_bam(file, peaks = x, genome = genome)
			cvg <- coverage(g)
			sum(cvg[x])
		}))

		se <- SummarizedExperiment(
			assays = list(counts = X),
			rowRanges = x,
			colData = data.frame(sample = names(filenames))
		)
		se
	}
) # get_counts

