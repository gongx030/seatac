#'
#' @export
#'
setMethod(
	'match_motifs',
	signature(
		x = 'Vplots',
		pwms = 'PWMatrixList'
	),
	function(x, pwms, genome){

		n_bins_per_window <- x@n_bins_per_window
		n_bins_per_block <- as.integer(x@block_size / x@bin_size)
		n_blocks_per_window <- n_bins_per_window - n_bins_per_block + 1

		motif_ix <- matchMotifs(homer_pwms, reduce(granges(x)), genome = genome, out = 'positions')

		classes <- names(motif_ix)
		tf_names <- gsub('(.+?)/.+', '\\1', classes)
		annotation <- unlist(motif_ix) %>%
			resize(fix = 'center', width = 1L)

		MS <- sparseMatrix(
			i = as.numeric(factor(names(annotation), classes)), 
			j = 1:length(annotation), 
			dims = c(length(classes), length(annotation))
		) # motifs ~ TF binding sites

		gr <- x %>% 
			granges() %>% 
			slidingWindows(x@block_size, x@bin_size) %>%
			unlist()

		gr <- gr %>% 
			resize(fix = 'center', width = x@bin_size)

		mm <- findOverlaps(annotation, gr, type = 'within') %>% 
		  as.matrix()

		SB <- sparseMatrix(i = mm[, 1], j = mm[, 2], dims = c(length(annotation), length(gr))) # TF binding site ~ block

		BM <- t(MS %*% SB)
		BM <- BM %>% as.array()
		dim(BM) <- c(length(x), n_blocks_per_window, length(classes))
    SummarizedExperiment::rowData(x)$motifs <- BM
		x
	}
)
	

