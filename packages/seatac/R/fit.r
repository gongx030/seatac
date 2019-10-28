#' fit.vae
#'
fit.vae <- function(model, gr, epochs = 1, batch_size = 256, learning_rate = 0.001){

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
  num_samples <- model$num_samples
	window_dim <- length(gr)

	x <- mcols(gr)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, model$feature_dim, model$input_dim, 1L))

	vae_loss <- function (y_true, y_pred) - (y_pred %>% tfd_log_prob(y_true))

  model$vae %>% compile(
		loss = vae_loss,
		optimizer = optimizer_adam(lr = learning_rate)
	)

	model$vae %>% fit(x, x, epochs = epochs, batch_size = batch_size)

	model

} # fit.vae


#' fit.cvae
#'
fit.cvae <- function(model, gr, epochs = 1, batch_size = 256, learning_rate = 0.001){

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
  num_samples <- model$num_samples
	window_dim <- length(gr)

	x <- mcols(gr)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, model$feature_dim, model$input_dim, 1L))

	sequence <- unlist(strsplit(as.character(mcols(gr)$sequence), '')) %>%
		factor(c('A', 'C', 'G', 'T')) %>%
		as.numeric() %>%
		matrix(nrow = window_dim, ncol = metadata(gr)$window_size, byrow = TRUE)

	sequence <- sequence - 1 	# convert to zero-based index

	vae_loss <- function (y_true, y_pred) - (y_pred %>% tfd_log_prob(y_true))

  model$vae %>% compile(
		loss = vae_loss,
		optimizer = optimizer_adam(lr = learning_rate)
	)

	model$vae %>% fit(list(x, sequence), x, epochs = epochs, batch_size = batch_size)

	model

} # fit.vae


#' fit.gmm_cvae
#'
fit.gmm_cvae <- function(model, gr, learning_rate = 0.001){

	batch_size <- metadata(gr)$training$batch_size
	steps_per_epoch <- metadata(gr)$training$steps_per_epoch
	epochs <- metadata(gr)$training$epochs

	flog.info(sprintf('batch size(batch_size): %d', batch_size))
	flog.info(sprintf('steps per epoch(steps_per_epoch): %d', steps_per_epoch))

	window_size <- metadata(gr)$window_size
	window_dim <- length(gr)

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))

	for (epoch in seq_len(epochs)) {

		total_loss <- 0
		total_loss_nll <- 0
		total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- gr$epoch == epoch & gr$step == s 	# window index for current batch

			x <- mcols(gr[b])$counts %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$feature_dim, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			sequence <- unlist(strsplit(as.character(mcols(gr[b])$sequence), '')) %>%
				factor(c('N', 'A', 'C', 'G', 'T')) %>%
				as.numeric() %>%
				matrix(nrow = batch_size, ncol = window_size, byrow = TRUE) %>%
				tf$cast(tf$float32)
			sequence <- sequence - 1

			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				g <- list(x, sequence) %>% model$encoder()

				posterior <- g[[1]]
				z_sequence <- g[[2]]

				posterior_sample <- posterior$sample()

				cond_posterior_sample <- posterior_sample + z_sequence

				likelihood <- cond_posterior_sample %>% model$decoder()

				nll <- -likelihood$log_prob(x)

				avg_nll <- tf$reduce_mean(nll)

				prior_model <- model$latent_prior_model(NULL)

				kl_div <- posterior$log_prob(posterior_sample) - prior_model$log_prob(posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)

				loss <- kl_div + avg_nll 
			})

			total_loss <- total_loss + loss
			total_loss_nll <- total_loss_nll + avg_nll
			total_loss_kl <- total_loss_kl + kl_div

			encoder_gradients <- tape$gradient(loss, model$encoder$variables)
			decoder_gradients <- tape$gradient(loss, model$decoder$variables)
			prior_gradients <- tape$gradient(loss, model$latent_prior_model$variables)

			optimizer$apply_gradients(
				purrr::transpose(list(encoder_gradients, model$encoder$variables)),
				global_step = tf$train$get_or_create_global_step()
			)

			optimizer$apply_gradients(
				purrr::transpose(list(decoder_gradients, model$decoder$variables)),
				global_step = tf$train$get_or_create_global_step()
			)

			optimizer$apply_gradients(
				purrr::transpose(list(prior_gradients, model$latent_prior_model$variables)),
				global_step = tf$train$get_or_create_global_step()
			)

		}

#		if (epoch == 1 || epoch %% 10 == 0){
#			gr2 <- model %>% predict(gr[sample.int(length(gr), 500)])
#			smoothScatter(gr2$latent, main = epoch)
#			points(prior_model$components_distribution$mean() %>% as.matrix(), bg = 1:4, pch = 21)
#		}

		flog.info(sprintf('training | epoch=%4.d/%4.d | nll=%7.1f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
	}
}
