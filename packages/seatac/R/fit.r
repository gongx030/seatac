#' fit.cvae
#'
fit.cvae <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)

	# compute the average vector at the fragment size dimension
	H <- metadata(gr)$n_bins_per_window *  metadata(gr)$n_intervals

	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$n_intervals), dims = c(H, metadata(gr)$n_bins_per_window))
	fragment_size <- mcols(gr)$counts %*% A %>% as.matrix()
	fragment_size <- Diagonal(x = 1 / rowSums(fragment_size)) %*% fragment_size %>% as.matrix()


	# compute the average vector at the position dimension
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_intervals, metadata(gr)$n_bins_per_window), dims = c(H, metadata(gr)$n_bins_per_window))
	position <- mcols(gr)$counts %*% A %>% as.matrix()

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			c_fragment_size <- fragment_size[b, ] %>% tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				g <- list(x, c_fragment_size) %>% model$encoder()
				posterior <- g[[1]]
				h_fragment_size <- g[[2]]

				posterior_sample <- posterior$sample()
				cond_posterior_sample <- tf$concat(list(posterior_sample, h_fragment_size), axis = 1L)

				likelihood  <- cond_posterior_sample %>% model$decoder()
				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	mcols(model$data)$fragment_size <- fragment_size
	model

} # fit.gmm_cvae

#' fit.vae
#'
fit.vae <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- Diagonal(x = 1 / rowSums(mcols(gr[b])$counts)) %*% mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae


#' fit.vae_baseline
#'
fit.vae_baseline <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae


#' fit.vae_20191216a
#'
fit.vae_20191216a <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			fragment_size <- x %>% tf$reduce_mean(axis = 2L) %>% tf$squeeze()
			position <- x %>% tf$reduce_mean(axis = 1L) %>% tf$squeeze()

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- list(fragment_size, position) %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae_20191216a

#' fit.vae_20191216b
#'
fit.vae_20191216b <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_vplot <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			fragment_size <- x %>% tf$reduce_mean(axis = 2L) %>% tf$squeeze()
			position <- x %>% tf$reduce_mean(axis = 1L) %>% tf$squeeze()

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- list(x, fragment_size, position) %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_vplot <- -likelihood[[1]]$log_prob(x)
				nll_fragment_size <- -likelihood[[2]]$log_prob(fragment_size)
				nll_position <- -likelihood[[3]]$log_prob(position)

				avg_nll_vplot <- tf$reduce_mean(nll_vplot)
				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_vplot + avg_nll_fragment_size + avg_nll_position
			})

			total_loss <- total_loss + loss
			total_loss_nll_vplot <- total_loss_nll_vplot + avg_nll_vplot
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(vplot)=%7.1f | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll_vplot, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae_20191216b


#' fit.vae_20191216c
#'
fit.vae_20191216c <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10, beta = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)


	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_vplot <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)
			x <- x + 1e-10

			fragment_size <- mcols(gr[b])$fragment_size %>% tf$cast(tf$float32) + 1e-5
			position <- mcols(gr[b])$position %>% tf$cast(tf$float32) + 1e-5

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_vplot <- -likelihood[[1]]$log_prob(x)
				nll_fragment_size <- -likelihood[[2]]$log_prob(fragment_size)
				nll_position <- -likelihood[[3]]$log_prob(position)

				avg_nll_vplot <- tf$reduce_mean(nll_vplot)
				avg_nll_fragment_size <- beta * tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- beta * tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_vplot + avg_nll_fragment_size + avg_nll_position
			})

			total_loss <- total_loss + loss
			total_loss_nll_vplot <- total_loss_nll_vplot + avg_nll_vplot
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(vplot)=%7.1f | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll_vplot, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae_20191216c


#' fit.vae_20191216d
#'
fit.vae_20191216d <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$latent_prior_model(NULL)


	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_vplot <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>% tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_vplot <- -likelihood[[1]]$log_prob(x)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_vplot <- tf$reduce_mean(nll_vplot)
				avg_nll_position <- tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_vplot + avg_nll_position
			})

			total_loss <- total_loss + loss
			total_loss_nll_vplot <- total_loss_nll_vplot + avg_nll_vplot
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(vplot)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll_vplot, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae_20191216d


#' fit.vae_20191216e
#'
fit.vae_20191216e <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10, beta = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	flog.info('computing fragment size vectors')
	H <- metadata(gr)$n_intervals * metadata(gr)$n_bins_per_window
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$n_intervals), dims = c(H, metadata(gr)$n_bins_per_window))
	fragment_size <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(x)$fragment_size <- (fragment_size - rowMins(fragment_size)) / (rowMaxs(fragment_size) - rowMins(fragment_size))

	flog.info('computing position vectors')
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_intervals, metadata(gr)$n_bins_per_window), dims = c(H, metadata(gr)$n_bins_per_window))
	position <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(x)$position <- (position - rowMins(position)) / (rowMaxs(position) - rowMins(position))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)
			x <- x + 1e-5

			fragment_size <- mcols(gr[b])$fragment_size %>% tf$cast(tf$float32) + 1e-5
			position <- mcols(gr[b])$position %>% tf$cast(tf$float32) + 1e-5

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_fragment_size <- beta * tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- beta * tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size + avg_nll_position
			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model
}


#' fit.vae_imputation
#'
fit.vae_imputation <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	flog.info('computing fragment size vectors')
	H <- metadata(gr)$n_intervals * metadata(gr)$n_bins_per_window
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_bins_per_window, each = metadata(gr)$n_intervals), dims = c(H, metadata(gr)$n_bins_per_window))
	fragment_size <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(gr)$fragment_size <- (fragment_size - rowMins(fragment_size)) / (rowMaxs(fragment_size) - rowMins(fragment_size))

	flog.info('computing position vectors')
	A <- sparseMatrix(i = 1:H, j = rep(1:metadata(gr)$n_intervals, metadata(gr)$n_bins_per_window), dims = c(H, metadata(gr)$n_bins_per_window))
	position <- mcols(gr)$counts %*% A %>% as.matrix()
	mcols(gr)$position <- (position - rowMins(position)) / (rowMaxs(position) - rowMins(position))

	prior_model <- model$latent_prior_model(NULL)

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			g <- (x == 0) %>% tf$dtypes$cast(tf$float32)	# dropout indicator

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				if (epoch == 1){
					posterior <- x %>% model$encoder()
				}else{
					z <- x %>% model$encoder()
					y <- z$mean() %>% model$decoder()	# decoded v-plot
					v <- y$mean() %>% model$imputer()
					xi <- x + v * g	# inputed v-plot
					posterior <- xi %>% model$encoder()
				}

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)
				avg_nll <- tf$reduce_mean(nll)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 

			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # vae_imputation
