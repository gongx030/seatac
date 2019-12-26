#' load_model
load_model <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	flog.info(sprintf('reading %s', model_file))
	x <- readRDS(model_file)

	model <- initialize_model(x)

	for (s in names(x))
		if (is.null(model[[s]]))
			model[[s]] <- x[[s]]
	model

} # load_model


initialize_model <- function(x, ...){

	UseMethod('initialize_model', x)

} # 

#' initialize_model.vae_baseline
#'
initialize_model.vae_baseline <- function(x){

	model <- build_model(
		type = x$type,
		input_dim = x$input_dim,
		feature_dim = x$feature_dim,
		latent_dim = x$latent_dim,
		num_samples = x$num_samples,
		window_size = x$window_size
	)

	vplot <- array(0, dim = c(1, model$feature_dim, model$input_dim, 1)) %>% 
		tf$cast(tf$float32)
	vplot %>% model$encoder()

	model$encoder$load_weights(x$encoder_weight_file)

	z <- matrix(0, nrow = 1, ncol = model$latent_dim) %>% 
		tf$cast(tf$float32)
	z %>% model$decoder()
	model$decoder$load_weights(x$decoder_weight_file)

	#model$prior$load_weights(x$latent_prior_weight_file)

	model

} # initialize_model.vae_baseline

