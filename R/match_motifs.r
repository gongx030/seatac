#'
#' @export
#'
setMethod(
	'match_motifs',
	signature(
		pwms = 'PWMatrixList',
		x = 'RangedSummarizedExperiment'
	),
	function(pwms, x, genome){

		motif_ix <- matchMotifs(pwms, reduce(resize(granges(x), fix = 'center', width = 50L)), genome = genome, out = 'position')

		classes <- names(motif_ix)
		tf_names <- sapply(pwms, function(p) p@name)
		annotation <- unlist(motif_ix) %>%
			resize(fix = 'center', width = 1L)

		MS <- sparseMatrix(
			i = as.numeric(factor(names(annotation), classes)), 
			j = 1:length(annotation), 
			dims = c(length(classes), length(annotation))
		) # motifs ~ TF binding sites

		mm <- findOverlaps(annotation, x, type = 'any') %>% 
		  as.matrix()

		SB <- sparseMatrix(i = mm[, 1], j = mm[, 2], dims = c(length(annotation), length(x))) # TF binding site ~ block

		BM <- t(MS %*% SB)
		BM <- BM %>% as('lgCMatrix')

		SummarizedExperiment(
			assays = list(motifMatches = BM),
			rowRanges = rowRanges(x),
			colData = data.frame(name = tf_names)
		)
	}
)
	

