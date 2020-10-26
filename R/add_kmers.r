#' add_kmers
#'
#' @export
#'
setMethod(
	'add_kmers',
	signature(
		x = 'Vplots',
		k = 'integer'
	),
	function(
		x,
		k,
		genome,
		...
	){
		gr <- granges(x) 
		res <- add_kmers_core(gr, k, genome)
		SummarizedExperiment::rowData(x)$kmers <- res$z
		class(x) <- 'VplotsKmers'
		x@kmers <- res$kmers
		x@k <- k
		x
	}
)


#' add_kmers_core
#'
add_kmers_core <- function(x, k, genome){

	window_size <- width(x)

	if (length(unique(window_size)) > 1)
		stop('the window size of input data must be equal')

	window_size <- window_size[1]

	nucleotides <- c('A', 'C', 'G', 'T')
	kmers <- do.call('paste0', do.call('expand.grid', rep(list(nucleotides), k)))

	start(x) <- start(x) - round(k / 2)
	end(x) <- end(x) + (k - round(k / 2) - 1)
	y <-  getSeq(genome, x) %>%
		as.character() %>%
		strsplit('')
	y <- do.call('rbind', y)

	invalid <- !y %in% nucleotides
	if (any(invalid)){
		y[invalid] <- sample(nucleotides, sum(invalid), replace = TRUE)
	}

	z <- lapply(1:window_size, function(start){
		do.call('paste0', as.data.frame(y[, start:(start + k - 1)])) %>%
			factor(kmers) %>%
			as.numeric()
	})
	z <- do.call('cbind', z)
	z <- z - 1 # change to zero-based for embedding
	class(z) <- 'integer'

	list(z = z, kmers = kmers)
} # add_kmers_core

