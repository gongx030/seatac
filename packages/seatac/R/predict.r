#' predict
#'
predict <- function(model, ...){

  if (model$input == 'fragment_size_position'){

		predict_vae_fragment_size_position(model, ...)

  }else
		stop('unknown model input')

} # predict


#' predict_vae_vplot_input
#'
predict_vae_vplot_input <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		g <- x %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict_vae_vplot_input



#' predict.cvae
#' 
predict.cvae <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	H <- metadata(gr)$n_bins_per_window *  metadata(gr)$n_intervals
	A <- sparseMatrix(
		i = 1:H, 
		j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$n_intervals), 
		dims = c(H, metadata(gr)$n_bins_per_window)
	)
	c_fragment_size <- mcols(gr)$counts %*% A %>% as.matrix()
	c_fragment_size <- Diagonal(x = 1 / rowSums(c_fragment_size)) %*% c_fragment_size %>% 
		as.matrix() %>%
		tf$cast(tf$float32)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z
	h_fragment_size <- matrix(NA, window_dim, model$hidden_dim)

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		g <- list(x, c_fragment_size[i, ]) %>% model$encoder()

		z[i, ] <- g[[1]]$mean() %>% as.matrix()
		h_fragment_size[i, ] <- g[[2]] %>% as.matrix()

	}

	mcols(gr)$latent <- z
	mcols(gr)$h_fragment_size <- h_fragment_size

	gr

} # predict.vae

#' predict.vae
#' 
predict.vae <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- Diagonal(x = 1 / rowSums(mcols(gr[i])$counts)) %*% mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		g <- x %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict.ae

#' predict.vae_baseline
#' 
predict.vae_baseline <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		g <- x %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict.vae_baseline


#' predict.vae_20191216a
#' 
predict.vae_20191216a <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		fragment_size <- x %>% tf$reduce_mean(axis = 2L) %>% tf$squeeze()
		position <- x %>% tf$reduce_mean(axis = 1L) %>% tf$squeeze()

		g <- list(fragment_size, position) %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict.vae_baseline


#' predict.vae_20191216b
#' 
predict.vae_20191216b <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		fragment_size <- x %>% tf$reduce_mean(axis = 2L) %>% tf$squeeze()
		position <- x %>% tf$reduce_mean(axis = 1L) %>% tf$squeeze()

		g <- list(x, fragment_size, position) %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict.vae_20191216b


#' predict.vae_20191216c
#' 
predict.vae_20191216c <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		g <- x %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict.vae_20191216c


#' predict.vae_20191216e
#' 
predict.vae_output_fragment_size_position <- function(model, gr, batch_size = 512){
	predict_vae_vplot_input(model = model, gr = gr, batch_size = batch_size)
}

#' predict.vae_imputation
#' 
predict.vae_imputation <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim)   # z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))

		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		posterior <- x %>% model$encoder()
		y <- posterior$mean() %>% model$decoder()
		v <- y$mean() %>% model$imputer()
		xi <- x + v # inputed v-plot
		posterior <- xi %>% model$encoder()

		z[i, ] <- posterior$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

}

#' predict.vae_imputation_gmm
#' 
predict.vae_imputation_gmm <- function(model, gr, batch_size = 512){
	predict_vae_vplot_input(model = model, gr = gr, batch_size = batch_size)
}

#' predict.vae_knn
#' 
predict.vae_knn <- function(model, gr, batch_size = 512){
	predict_vae_vplot_input(model = model, gr = gr, batch_size = batch_size)
}

#' predict.vae_output_fragment_size_position_with_imputation
#' 
predict.vae_output_fragment_size_position_with_imputation <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim)   # z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))

		i <- starts[b]:ends[b]

		x <- mcols(gr[i])$counts  %>%
			as.matrix() %>%
			array_reshape(c(length(i), model$feature_dim, model$input_dim, 1L)) %>%
			tf$cast(tf$float32)

		zero_entry <- (x == 0) %>% tf$cast(tf$float32)

		posterior <- x %>% model$encoder()
		y <- posterior$mean() %>% model$decoder()
		v <- list(y[[1]]$sample(), y[[2]]$sample()) %>% model$imputer()
		xi <- x + 0.1 * v * zero_entry  # inputed v-plot
		posterior <- xi %>% model$encoder()

		z[i, ] <- posterior$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # vae_output_fragment_size_position_with_imputation


#' predict_vae_fragment_size_position
#'
predict_vae_fragment_size_position <- function(model, gr, batch_size = 512){

	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size
	num_samples <- metadata(gr)$num_samples
	n_bins_per_window <- metadata(gr)$window_size / metadata(gr)$bin_size
	bin_size <- metadata(gr)$bin_size

	flog.info(sprintf('# input peaks: %d', window_dim))
	flog.info(sprintf('batch size(batch_size): %d', batch_size))

	starts <- seq(1, window_dim, by = batch_size)
	ends <- starts + batch_size - 1
	ends[ends > window_dim] <- window_dim
	n_batch <- length(starts)

	z <- matrix(NA, window_dim, model$latent_dim) 	# z

	for (b in 1:n_batch){

		if (b %% 10 == 0)
			flog.info(sprintf('predicting | batch=%4.d/%4.d', b, n_batch))
		
		i <- starts[b]:ends[b]

		fragment_size <- mcols(gr[i])$fragment_size %>% 
			as.matrix() %>%
			tf$cast(tf$float32)

		position <- mcols(gr[i])$position %>% 
			as.matrix() %>% 
			tf$cast(tf$float32)

		g <- list(fragment_size, position) %>% model$encoder()
		z[i, ] <- g$mean() %>% as.matrix()

	}

	mcols(gr)$latent <- z

	gr

} # predict_vae_fragment_size_position


