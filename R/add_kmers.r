add_kmers <- function(x, kmer, kmer_step){

	nucleotides <- c('A', 'C', 'G', 'T')
	kmers <- expand.grid(lapply(1:kmer, function(i) nucleotides), stringsAsFactors = FALSE)
	kmers <- Reduce('paste0', kmers)

	sequence <- as.character(mcols(x)$sequence)
	starts <- seq(1, metadata(x)$window_size, by = kmer_step)

	sequence <- unlist(lapply(starts, function(start) substr(sequence, start, start + kmer - 1)))
	invalid <- !sequence %in% kmers

	# string distance between invalid kmers (containing Ns) with all valid kmers
	D <- adist(sequence[invalid], kmers)
	sequence[invalid] <- kmers[max.col(-D)]
	sequence <- sequence %>% 
		factor(kmers) %>% 
		as.numeric() %>%
		matrix(nrow = length(x), ncol = length(starts))

	mcols(x)$kmers <- sequence
	metadata(x)$kmers <- kmers
	metadata(x)$kmer_width <- kmer
	metadata(x)$kmer_step <- kmer_step
	metadata(x)$kmer_starts <- starts
	x

} # add_kmers
