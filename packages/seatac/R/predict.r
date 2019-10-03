#' predict.vae
#' 
predict.vae <- function(model, gr, chunk_size = 2^17, batch_size = 256){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	starts <- seq(1, window_dim, by = chunk_size)
	ends <- starts + chunk_size - 1
	ends[ends > window_dim] <- window_dim
	n_chunk <- length(starts)

	# latent representation of each unique window
	Z <- matrix(NA, window_dim, model$latent_dim)

	for (b in 1:n_chunk){

		flog.info(sprintf('%d/%d', b, n_chunk))
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L))

		Z[i, ]  <- model$latent %>% predict(x, batch_size = batch_size, verbose = 1)
	}

	mcols(gr)$latent <- Z

	gr

} # predict.vae


#' predict.cvae
#' 
predict.cvae <- function(model, gr, chunk_size = 2^17, batch_size = 256){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))

	starts <- seq(1, window_dim, by = chunk_size)
	ends <- starts + chunk_size - 1
	ends[ends > window_dim] <- window_dim
	n_chunk <- length(starts)

	# latent representation of each unique window
	z <- matrix(NA, window_dim, model$latent_dim) 	# z
	nfr <- matrix(NA, window_dim, model$input_dim)
	mono_nucleosome <- matrix(NA, window_dim, model$input_dim)

	for (b in 1:n_chunk){

		flog.info(sprintf('%d/%d', b, n_chunk))
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L))

		sequence <- unlist(strsplit(as.character(mcols(gr[i])$sequence), '')) %>%
			factor(c('A', 'C', 'G', 'T')) %>%
			as.numeric() %>%
			matrix(nrow = length(i), ncol = window_size, byrow = TRUE)

		sequence <- sequence - 1

		h <- model$latent %>% predict(list(x, sequence), batch_size = batch_size, verbose = 1)
		z[i, ] <- h[[1]]
		nfr[i, ] <- h[[2]][, , 1]
		mono_nucleosome[i, ] <- h[[3]][, , 1]
	}

	mcols(gr)$latent <- z
	mcols(gr)$nfr <- nfr
	mcols(gr)$mono_nucleosome <- mono_nucleosome

	gr

} # predict.cvae


#' predict.gmm_cvae
#' 
predict.gmm_cvae <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	# latent representation of each unique window
	z <- matrix(NA, window_dim, model$latent_dim) 	# z
	nfr <- matrix(NA, window_dim, model$input_dim)
	mono_nucleosome <- matrix(NA, window_dim, model$input_dim)

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		sequence <- unlist(strsplit(as.character(mcols(gr[i])$sequence), '')) %>%
			factor(c('N', 'A', 'C', 'G', 'T')) %>%
			as.numeric() %>%
			matrix(nrow = length(i), ncol = window_size, byrow = TRUE)

		sequence <- sequence - 1

		z[i, ] <- model$encoder(list(x, sequence))[[1]]$mean() %>% as.matrix()
	}

	mcols(gr)$latent <- z

	gr

} # predict.gmm_cvae
