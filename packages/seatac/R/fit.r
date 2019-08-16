#' fit.vae
#'
fit.vae <- function(model, gr, epochs = 1, batch_size = 256, learning_rate = 0.001){

	optimizer <- tf$train$AdamOptimizer(learning_rate)
	flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))
  num_samples <- model$num_samples
	window_dim <- length(gr)

	x <- mcols(gr)$counts %>%
		as.matrix() %>%
		array_reshape(c(window_dim, model$input_dim, model$feature_dim, 1L)) 

	y <- mcols(gr)$coverage %>%
		array_reshape(c(window_dim, model$input_dim, 1L))

	label <- as.numeric(mcols(gr)$label)

	vae_loss <- function (y_true, y_pred) - (y_pred %>% tfd_log_prob(y_true))

  model$vae %>% compile(
		loss = vae_loss,
		optimizer = optimizer_adam(lr = learning_rate)
	)

	model$vae %>% fit(list(x, y), list(x, y, label), epochs = epochs, batch_size = batch_size)

	model

} # fit_vae
