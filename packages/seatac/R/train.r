fit_vae <- function(model, train_dataset, epochs = 1, steps_per_epoch = 10){

	optimizer <- tf$train$AdamOptimizer(0.001)
	flog.trace('optimizer: Adam(learning_rate=0.001)')
	flog.trace(sprintf('steps_per_epoch: %d', steps_per_epoch))
	flog.trace(sprintf('trainable prior: %s', model$trainable_prior))

  for (epoch in seq_len(epochs)) {

    iter <- make_iterator_one_shot(train_dataset)

    total_loss <- 0
    total_loss_nll <- 0
    total_loss_kl <- 0

		for (s in 1:steps_per_epoch){

      x <- iterator_get_next(iter)

      with(tf$GradientTape(persistent = TRUE) %as% tape, {

        approx_posterior <- model$encoder(x)
        approx_posterior_sample <- approx_posterior$sample()
        decoder_likelihood <- model$decoder(approx_posterior_sample)

        nll <- -decoder_likelihood$log_prob(x)
        avg_nll <- tf$reduce_mean(nll)

				if (model$trainable_prior)
					model$latent_prior <- model$latent_prior_model(NULL)

				kl_div <- approx_posterior$log_prob(approx_posterior_sample) - model$latent_prior$log_prob(approx_posterior_sample)
				kl_div <- tf$reduce_mean(kl_div)
        loss <- kl_div + avg_nll
      })

      total_loss <- total_loss + loss
      total_loss_nll <- total_loss_nll + avg_nll
      total_loss_kl <- total_loss_kl + kl_div

      encoder_gradients <- tape$gradient(loss, model$encoder$variables)
      decoder_gradients <- tape$gradient(loss, model$decoder$variables)

      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients, model$encoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients, model$decoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

			if (model$trainable_prior){
      	prior_gradients <- tape$gradient(loss, model$latent_prior_model$variables)
				optimizer$apply_gradients(
					purrr::transpose(list(prior_gradients, model$latent_prior_model$variables)),
					global_step = tf$train$get_or_create_global_step()
				)
			}
    }

    flog.info(sprintf('epoch=%d/%d | negative likelihood=%.3f | kl=%.3f | total=%.3f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
  }

} # fit_vae
