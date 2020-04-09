setGeneric('add_gc_content', function(x, genome, ...) standardGeneric('add_gc_content'))

#' add_gc_content
#'
#' @export
#'
setMethod(
	'add_gc_content',
	signature(
		x = 'GRanges',
		genome = 'BSgenome'
	),
	function(x, genome, ...){
		mcols(x)$gc_content <- getSeq(genome, x) %>%
			letterFrequency(letters = 'CG', as.prob = TRUE)
		x
	}
)

