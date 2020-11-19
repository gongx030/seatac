#' Compute deviations
#' 
#' Compute the V-plot deviations
#'
#' @param x a Vplots object 
#' @param annotation a GRangesList object of motif binding sites
#' @param model a pretrained VaeModel
#' @param batch_size_window batch size for processing windows (default: 32L)
#' @param batch_size_block batch size for running VAE prediction (default: 128L)
#' @param background number of background V-plots
#' @param permutation number of permutations
#'
#' @return a SummarizedVplots object
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'compute_deviations',
	signature(
		x = 'RangedSummarizedExperiment',
		annotation = 'RangedSummarizedExperiment'
	), 
	function(
		x, 
		annotation,
		permutation = 100L
	){

		BM <- assays(annotation)$motifMatches
		BC <- assays(x)$counts

		BCE <- rowSums(BC) %*% t(colSums(BC)) / sum(BC)
		Y <- (t(BM) %*% BC - t(BM) %*% BCE) / (t(BM) %*% BCE) # raw deviation
		Y <- Y %>% as.matrix() %>% tf$cast(tf$float32) %>% tf$expand_dims(0L)

		Y_perm <- lapply(1:permutation, function(i){
			j <- sample(1:length(x))
			y <- (t(BM[j, ]) %*% BC - t(BM[j, ]) %*% BCE) / (t(BM) %*% BCE) # raw deviation
			y %>% as.matrix() %>% tf$cast(tf$float32) %>% tf$expand_dims(0L)
		})

		Y_perm <- tf$concat(Y_perm, axis = 0L)

		Y_mean <- Y_perm %>% tf$reduce_mean(0L, keepdims = TRUE)	# mean nucleosome
		Y_std <- Y_perm %>% tf$math$reduce_std(0L, keepdims = TRUE)	# std nucleosome

		Z <- (Y - Y_mean) / Y_std
		Z <- Z %>% tf$squeeze(0L)
		Z <- Z %>% as.matrix()
		rownames(Z) <- colnames(annotation)
		browser()


		browser()

	}
) # compute_deviations
