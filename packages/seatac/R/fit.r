fit.vae <- function(model, gr, epochs = 1, steps_per_epoch = 10, batch_size = 256, beta = 1, gamma = 10, learning_rate = 0.001){

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
	flog.info(sprintf('trainable prior: %s', model$trainable_prior))
  num_samples <- model$num_samples


  for (epoch in seq_len(epochs)) {

		# randomly sampling windows
		G <- matrix(sample.int(length(gr), steps_per_epoch * batch_size), steps_per_epoch, batch_size)

    total_loss <- 0
    total_loss_nll_x <- 0
    total_loss_nll_y <- 0
    total_loss_nll_label <- 0
    total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

			b <- G[s, ]	# window index for current batch
			x <- mcols(gr)$counts[b, ] %>%
				as.matrix() %>%
				array_reshape(c(batch_size, model$input_dim, model$feature_dim, 1L)) %>%
				tf$cast(tf$float32)	# need to convert to tensors

#			y <- mcols(gr)$coverage[b, , drop = FALSE] %>%
#				array_reshape(c(batch_size, model$input_dim, 1L)) %>%
#				tf$cast(tf$float32)

			label <- as.numeric(mcols(gr)$label[b])
			
			y <- mcols(gr)$coverage[b, , drop = FALSE]; y[, 1:10] <- label
			y <- y %>%
				array_reshape(c(batch_size, model$input_dim, 1L)) %>%
				tf$cast(tf$float32)

			# determining the weight for each window
			w <- windows_weight(mcols(gr)$num_reads[b])

      with(tf$GradientTape(persistent = TRUE) %as% tape, {

        posterior_x <- x %>% model$encoder$vplot()
        posterior_sample_x <- posterior_x$sample()

        posterior_y <- y %>% model$encoder$coverage()
        posterior_sample_y <- posterior_y$sample()

				label_dist <- list(vplot = x, coverage = y) %>% model$classifier()
				label_pred <- label_dist$sample()

				posterior_sample <- tf$concat(list(posterior_sample_x, posterior_sample_y, label_pred %>% tf$cast(tf$float32)), 1L)

       	likelihood_x <- posterior_sample %>% model$decoder$vplot()
       	likelihood_y <- posterior_sample %>% model$decoder$coverage()

				nll_x <- -w * likelihood_x$log_prob(x)
				nll_y <- -w * likelihood_y$log_prob(y)
				nll_label <- - label_dist$log_prob(label)

				tab <- table(label, label_pred %>% as.matrix())
				print((tab[1, 1] + tab[2, 2]) / sum(tab))

        avg_nll_x <- tf$reduce_mean(nll_x)
        avg_nll_y <- tf$reduce_mean(nll_y)
        avg_nll_label <- gamma * tf$reduce_mean(nll_label)

				prior_model <- model$latent_prior_model(NULL)

				if (model$prior == 'gmm'){
					kl_div <- w * (posterior_x$log_prob(posterior_sample_x) + posterior_y$log_prob(posterior_sample_y)  - prior_model$log_prob(posterior_sample))
				}

				kl_div <- beta * tf$reduce_mean(kl_div)

        loss <- 0 * kl_div + 0 * avg_nll_x + 0 * avg_nll_y + avg_nll_label
      })

      total_loss <- total_loss + loss
			total_loss_nll_x <- total_loss_nll_x + avg_nll_x
      total_loss_nll_y <- total_loss_nll_y + avg_nll_y
      total_loss_nll_label <- total_loss_nll_label + avg_nll_label
      total_loss_kl <- total_loss_kl + kl_div

      encoder_gradients_x <- tape$gradient(loss, model$encoder$vplot$variables)
      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients_x, model$encoder$vplot$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      encoder_gradients_y <- tape$gradient(loss, model$encoder$coverage$variables)
      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients_y, model$encoder$coverage$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      decoder_gradients_x <- tape$gradient(loss, model$decoder$vplot$variables)
      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients_x, model$decoder$vplot$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      decoder_gradients_y <- tape$gradient(loss, model$decoder$coverage$variables)
      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients_y, model$decoder$coverage$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

			classifier_gradients <- tape$gradient(loss, model$classifier$variables)
      optimizer$apply_gradients(
        purrr::transpose(list(classifier_gradients, model$classifier$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

     	prior_gradients <- tape$gradient(loss, model$latent_prior_model$variables)
			optimizer$apply_gradients(
				purrr::transpose(list(prior_gradients, model$latent_prior_model$variables)),
				global_step = tf$train$get_or_create_global_step()
			)
    }

    flog.info(sprintf('training | epoch=%4.d/%4.d | nll(vplot)=%7.1f | nll(coverage)=%7.1f | nll(classifier)=%7.5f | kl=%7.1f | total=%7.1f', epoch, epochs, total_loss_nll_x, total_loss_nll_y, total_loss_nll_label, total_loss_kl, total_loss))
  }

} # fit_vae
