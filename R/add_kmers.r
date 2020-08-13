setMethod(
	'add_kmers',
	signature(
		x = 'Vplots',
		k = 'integer',
		genome = 'BSgenome'
	),
	function(
		x,
		k,
		genome,
		...
	){

		nucleotides <- c('A', 'C', 'G', 'T')
		kmers <- do.call('paste0', do.call('expand.grid', rep(list(nucleotides), k)))

		gr <- granges(x) 
		start(gr) <- start(gr) - round(k / 2)
		end(gr) <- end(gr) + (k - round(k / 2) - 1)
		y <-  getSeq(genome, gr) %>%
			as.character() %>%
			strsplit('')
		y <- do.call('rbind', y)

		invalid <- !y %in% nucleotides
		if (any(invalid)){
			y[invalid] <- sample(nucleotides, sum(invalid), replace = TRUE)
		}

		z <- lapply(1:x@window_size, function(start){
			do.call('paste0', as.data.frame(y[, start:(start + k - 1)])) %>%
				factor(kmers) %>%
				as.numeric()
		})
		z <- do.call('cbind', z)
		
		mcols(x)$kmers <- z - 1 # change to zero-based for embedding
		class(x) <- 'VplotsKmers'
		x@kmers <- kmers
		x@k <- k
		
		x
	}
)
