#' fit
#'
fit <- function(model, ...){

	if (model$input == 'fragment_size_position' && model$output == 'fragment_size_position' && !model$imputation && is.null(model$regularizer)){

		fit_fragment_size_position_input_fragment_size_position_output(model, ...)

	}else if (model$input == 'fragment_size_position' && model$output == 'fragment_size_position' && model$imputation && is.null(model$regularizer)){

		fit_fragment_size_position_input_fragment_size_position_output_with_imputation(model, ...)

	}else if (model$input == 'fragment_size' && model$output == 'vplot_parametric' && !model$imputation && is.null(model$regularizer)){
		fit_fragment_size_input_vplot_parametric_output(model, ...)
	}else if (model$input == 'fragment_size' && model$output == 'fragment_size' && !model$imputation && is.null(model$regularizer)){
		fit_fragment_size(model, ...)
	}else if (model$input == 'fragment_size_position' && model$output == 'fragment_size_position' && model$regularizer$name == 'kmeans'){
		fit_fragment_size_position_kmeans(model, ...)
	}else
		stop('unknown model input/output')

} # fit

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

	prior_model <- model$prior(NULL)

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

	prior_model <- model$prior(NULL)

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

	prior_model <- model$prior(NULL)

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

		tf$reduce_sum(x, axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'input data')
		tf$reduce_sum(likelihood$mean(), axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'decoded data')

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

	prior_model <- model$prior(NULL)

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

	prior_model <- model$prior(NULL)

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

	prior_model <- model$prior(NULL)


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

	prior_model <- model$prior(NULL)


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


#' fit.vae_output_fragment_size_position
#'
fit.vae_output_fragment_size_position <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$prior(NULL)

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

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>% 
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>% 
				as.matrix() %>% 
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- tf$reduce_mean(nll_position)

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

} # vae_output_fragment_size_position


#' fit.vae_imputation
#'
fit.vae_imputation <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10, zero_weight = 1e-3){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))
	flog.info(sprintf('setting weight for non-zero entries to 1 and zero entries to %.3e', zero_weight))


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

			zero_entry <- (x == 0) %>% tf$cast(tf$float32)

			# it is likely that too sparse will cause NA during the training
			weight <- (x > 0) %>% tf$dtypes$cast(tf$float32) + zero_weight # dropout indicator

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				if (epoch < 10){
					posterior <- x %>% model$encoder()
					posterior_sample <- posterior$sample()
					likelihood  <- posterior_sample %>% model$decoder()
					nll <- -likelihood$log_prob(x)
				}else{
					status <- 'imputation'
					z <- x %>% model$encoder()
					y <- z$sample() %>% model$decoder()	
					v <- y$sample() %>% model$imputer()
					xi <- x + v * zero_entry	# inputed v-plot
					posterior <- xi %>% model$encoder()
					posterior_sample <- posterior$sample()
					likelihood  <- posterior_sample %>% model$decoder()
					nll <- -likelihood$log_prob(x)
				}

				avg_nll <- tf$reduce_mean(nll)

				prior_model <- model$prior(NULL)
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

		tf$reduce_sum(x, axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'input data')
		tf$reduce_sum(likelihood$mean(), axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'decoded data')

		if (status == 'imputation'){
			tf$reduce_sum(xi, axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'imputed data')
		}

		flog.info(sprintf('%s | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', status, epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # vae_imputation


#' fit.vae_imputation_gmm
#'
fit.vae_imputation_gmm <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10, zero_weight = 1e-3){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))
	flog.info(sprintf('setting weight for non-zero entries to 1 and zero entries to %.3e', zero_weight))

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

			zero_entry <- (x == 0) %>% tf$dtypes$cast(tf$float32)	# dropout indicator

			# it is likely that too sparse will cause NA during the training
			weight <- (x > 0) %>% tf$dtypes$cast(tf$float32) + zero_weight	# dropout indicator

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				if (epoch == 1){
					status <- 'warming up'
					posterior <- x %>% model$encoder()
				}else{
					status <- 'imputation'
					z <- x %>% model$encoder()
					y <- z$mean() %>% model$decoder()	# decoded v-plot
					v <- y$mean() %>% model$imputer()
					xi <- x + v * zero_entry	# inputed v-plot
					posterior <- xi %>% model$encoder()
				}

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x) * weight 
				avg_nll <- tf$reduce_mean(nll)

				prior_model <- model$prior(NULL)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 

			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)
			prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))

		}

		flog.info(sprintf('%s | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', status, epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # vae_imputation_gmm


#' fit.vae_knn
#'
fit.vae_knn <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))
	
	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			c <- mcols(gr[b])$counts 

			x <- c %>% 
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- x %>% model$encoder()

				nn <- posterior$sample() %>% 
					as.matrix() %>% 
					knn.index(k = model$k)
				A <- sparseMatrix(i = rep(1:batch_size, model$k), j = c(nn), dims = c(batch_size, batch_size))

				xi <- A %*% c %>%
					as.matrix() %>%
					array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
					tf$cast(tf$float32) 
					
				posterior <- xi %>% model$encoder()
				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(xi)
				avg_nll <- tf$reduce_mean(nll)

				prior_model <- model$prior(NULL)

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

} # vae_knn



#' fit.vae_output_fragment_size_position_with_imputation
#'
fit.vae_output_fragment_size_position_with_imputation <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	prior_model <- model$prior(NULL)

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

			zero_entry <- (x == 0) %>% tf$cast(tf$float32)

			fragment_size <- mcols(gr[b])$fragment_size %>% tf$cast(tf$float32) + 1e-5
			position <- mcols(gr[b])$position %>% tf$cast(tf$float32) + 1e-5

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				z <- x %>% model$encoder()
				y <- z$sample() %>% model$decoder()	
				v <- list(y[[1]]$sample(), y[[2]]$sample()) %>% model$imputer()
				xi <- x + 0.1 * v * zero_entry	# inputed v-plot
				posterior <- xi %>% model$encoder()
				posterior_sample <- posterior$sample()
				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- tf$reduce_mean(nll_position)

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

		tf$reduce_sum(x, axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'input data')
		tf$reduce_sum(xi, axis = 0L) %>% tf$squeeze() %>% as.matrix() %>% t() %>% image(main = 'imputed data')

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit.vae_output_fragment_size_position_with_imputation


#' fit_fragment_size_position_input_fragment_size_position_output
#'
fit_fragment_size_position_input_fragment_size_position_output <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- list(fragment_size, position) %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)

				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size + avg_nll_position

			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

			print(model$decoder$trainable_variables)

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))
	}

	model$data <- gr
	model

} # fit_fragment_size_position_input_fragment_size_position_output


#' fit_fragment_size_position_input_fragment_size_position_output_with_imputation
#'
fit_fragment_size_position_input_fragment_size_position_output_with_imputation <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			fragment_size_zero <- (fragment_size == 0) %>% tf$cast(tf$float32)
			position_zero <- (position == 0) %>% tf$cast(tf$float32)

			if (epoch == 1){
				fragment_size_imputed <- fragment_size
				position_imputed <- position
			}else{
				z <- list(fragment_size, position) %>% model$encoder()
				h <- z$mean() %>% model$decoder()
				fragment_size_decoded <- h$fragment_size$mean()
			}

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

					
				posterior <- list(fragment_size_imputed, position_imputed) %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size_imputed) * (1 - fragment_size_zero)
				nll_position <- -likelihood[[2]]$log_prob(position_imputed) * (1 - position_zero)

				avg_nll_fragment_size <- tf$reduce_sum(nll_fragment_size)
				avg_nll_position <- tf$reduce_sum(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size + avg_nll_position
			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss))

		browser()
	}

	model$data <- gr
	model

} # fit_fragment_size_position_input_fragment_size_position_output_with_imputation


#' fit_fragment_size_input_vplot_parametric_output
#'
fit_fragment_size_input_vplot_parametric_output <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- fragment_size %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood$log_prob(fragment_size)
				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size 
			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

		}

		model$decoder$trainable_variables[[3]][1, ] %>% as.matrix() %>% plot(main = epoch)
		model$decoder$trainable_variables[[3]][2, ] %>% as.matrix() %>% plot(main = epoch)

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, total_loss_kl, total_loss))

	}

	model$data <- gr
	model

} # fit_fragment_size_input_vplot_parametric_output


#' fit_fragment_size
#'
fit_fragment_size <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				posterior <- fragment_size %>% model$encoder()

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood$log_prob(fragment_size)

				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)

				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size
			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

		}

#		gr2 <- model %>% predict(gr)
#		smoothScatter(gr2$latent, main = epoch)
#		model$prior(NULL)$components_distribution$mean() %>% as.matrix() %>% points(pch =21, bg = 'red')

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f kl=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, total_loss_kl, total_loss))
	}

	model$data <- gr
	model
} # fit_fragment_size



#' fit_fragment_size_position_kmeans
#'
fit_fragment_size_position_kmeans <- function(model, gr, learning_rate = 0.001, batch_size = 128, epochs = 50, steps_per_epoch = 10){

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)	# number of samples

	optimizer <- tf$keras$optimizers$Adam(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	if (length(model$prior$trainable_variables) > 0)
		trainable_prior <- TRUE
	else
		trainable_prior <- FALSE

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll_fragment_size <- 0
		total_loss_nll_position <- 0
		total_loss_kl <- 0
		total_loss_reg <- 0

		for (s in 1:steps_per_epoch){

			b <- sample.int(window_dim, batch_size)

			fragment_size <- mcols(gr[b])$fragment_size %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			position <- mcols(gr[b])$position %>% 
				as.matrix() %>%
				tf$cast(tf$float32)

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				h <- list(fragment_size, position) %>% model$encoder()

				posterior <- h[[1]]	# the first is always the latent distributioj

				posterior_sample <- posterior$sample()

				likelihood  <- posterior_sample %>% model$decoder()

				nll_fragment_size <- -likelihood[[1]]$log_prob(fragment_size)
				nll_position <- -likelihood[[2]]$log_prob(position)

				avg_nll_fragment_size <- tf$reduce_mean(nll_fragment_size)
				avg_nll_position <- tf$reduce_mean(nll_position)

				kl_div <- posterior$log_prob(posterior_sample) - model$prior(NULL)$log_prob(posterior_sample)

				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll_fragment_size + avg_nll_position

				if (!is.null(model$regularizer)){

					if (model$regularizer$name == 'kmeans'){
						z <- posterior$mean()
						reg_loss <- list(z, h[['mixture']]) %>% model$regularizer()
					}

					loss <- loss + reg_loss

				}

			})

			total_loss <- total_loss + loss
			total_loss_nll_fragment_size <- total_loss_nll_fragment_size + avg_nll_fragment_size
			total_loss_nll_position <- total_loss_nll_position + avg_nll_position
			total_loss_kl <- total_loss_kl + kl_div

			if (!is.null(model$regularizer)){
				total_loss_reg <- total_loss_reg + reg_loss
			}

			encoder_gradients <- tape$gradient(loss, model$encoder$trainable_variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$trainable_variables)

			if (trainable_prior){
				prior_gradients <- tape$gradient(loss, model$prior$trainable_variables)
			}

			if (!is.null(model$regularizer)){
				regularizer_gradients <- tape$gradient(loss, model$regularizer$trainable_variables)
			}

			optimizer$apply_gradients(purrr::transpose(list(encoder_gradients, model$encoder$trainable_variables)))
			optimizer$apply_gradients(purrr::transpose(list(decoder_gradients, model$decoder$trainable_variables)))

			if (trainable_prior){
				optimizer$apply_gradients(purrr::transpose(list(prior_gradients, model$prior$trainable_variables)))
			}

			if (!is.null(model$regularizer)){
				optimizer$apply_gradients(purrr::transpose(list(regularizer_gradients, model$regularizer$trainable_variables)))
			}

		}

		gr2 <- model %>% predict(gr)
		smoothScatter(gr2$latent, main = epoch)
#		model$prior(NULL)$components_distribution$mean() %>% as.matrix() %>% points(pch =21, bg = 'red')
#		model$0(NULL)$components_distribution$mean() %>% as.matrix() %>% points(pch =21, bg = 'red')
		model$regularizer$trainable_variables[[1]][, ] %>% as.matrix() %>% points(pch = 21, bg = 'red')

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll(fragment_size)=%7.1f | nll(position)=%7.1f | kl=%7.1f | reg=%7.1f | total=%7.1f', epoch, epochs, avg_nll_fragment_size, avg_nll_position, total_loss_kl, total_loss_reg, total_loss))
	}

	model$data <- gr
	model

} # fit_fragment_size_position_kmeans
