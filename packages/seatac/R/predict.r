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

		g <- model$encoder(list(x, sequence))
		z[i, ] <- g[[1]]$mean() %>% as.matrix()
		h <- g[[2]] %>% as.matrix()
		zh <- model$decoder(z[i, ] + h)$mean()

		nfr[i, ] <- zh[, which(metadata(gr)$nfr), , ] %>% 
			tf$reduce_mean(axis = 1L) %>%
			tf$squeeze() %>%
			as.matrix()

		mono_nucleosome[i, ] <- zh[, which(metadata(gr)$mono_nucleosome), , ] %>% 
			tf$reduce_mean(axis = 1L) %>%
			tf$squeeze() %>%
			as.matrix()
		
	}

	mcols(gr)$latent <- z
	mcols(gr)$nfr <- nfr
	mcols(gr)$mono_nucleosome <- mono_nucleosome

	n_bins <- 10
	core_bins <- 4
	weight <- exp(-((-n_bins:n_bins) / core_bins)^2 / 2)

	R <- log(mono_nucleosome * mean(nfr) / mean(mono_nucleosome) + 1e-5) - log(nfr + 1e-5)
	R <- do.call('cbind', lapply(1:n_bins_per_window, function(j){
		js <- (j - n_bins):(j + n_bins)
		valid <- js > 0 & js < n_bins_per_window
		rowSums(R[, js[valid]] %*% diag(weight[valid]))
	}))
	Z <- (R - mean(R)) / sd(c(R))
	mcols(gr)$seatac_nucleosome_ratio <- R
	mcols(gr)$seatac_nucleosome_score <- Z
	mcols(gr)$seatac_nucleosome_log10pvalue <- -log10(1 - pnorm(Z) + 1e-20)

	gr

} # predict.gmm_cvae
