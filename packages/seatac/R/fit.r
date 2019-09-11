#' fit.vae
#'
fit.vae <- function(model, gr, epochs = 1, batch_size = 256, learning_rate = 0.001){

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
  num_samples <- model$num_samples
	window_dim <- length(gr)
	window_size <- metadata(gr)$window_size

	S <- matrix(1:window_dim, window_dim / num_samples, num_samples)

	x_list <- lapply(1:num_samples, function(i){
		mcols(gr[S[, i]])$counts %>%
			as.matrix() %>%
			array_reshape(c(window_dim / num_samples, model$feature_dim, model$input_dim, 1L)) 
	})

	vae_loss <- function (y_true, y_pred) - (y_pred %>% tfd_log_prob(y_true))

  model$vae %>% compile(
		loss = vae_loss,
		optimizer = optimizer_adam(lr = learning_rate)
	)

	model$vae %>% fit(x_list, x_list, epochs = epochs, batch_size = batch_size)

	model

} # fit.vae
