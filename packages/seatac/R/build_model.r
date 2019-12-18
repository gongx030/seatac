#' build_model
#'
build_model <- function(type, input_dim, feature_dim, latent_dim, num_samples, window_size, ...){

	if (type == 'cvae'){

	}else if(type == 'vae'){

	}else if(type == 'vae_baseline'){

		structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_baseline_dense(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size
		), class = type)

	}else if(type == 'vae_20191216a'){
	}else if(type == 'vae_20191216b'){
	}else if(type == 'vae_20191216c'){
	}else if(type == 'vae_20191216d'){
	}else if(type == 'vae_output_fragment_size_position'){

		 structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size
		), class = c('vae_output_fragment_size_position'))

	}else if(type == 'vae_imputation'){

		#' This model attemps to impute the zero entries in the V-plot during the training
		structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_baseline_conv2(input_dim, feature_dim),
			imputer = imputer_model_baseline(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size
		), class = type)

	}else if(type == 'vae_imputation_gmm'){

		#' This model attemps to impute the zero entries in the V-plot during the training
		structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_baseline2(input_dim, feature_dim),
			imputer = imputer_model_baseline(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size
		), class = type)

	}else if(type == 'vae_knn'){

		#' This model attemps to impute the zero entries in the V-plot during the training
		structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_baseline(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			k = 20
		), class = type)

	}else if(type == 'vae_output_fragment_size_position_with_imputation'){

		 structure(list(
			encoder = encoder_model_vae_baseline(latent_dim, window_size),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			imputer = imputer_model_fragment_size_position(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size
		), class = type)

	}else if(type == 'vae_fragment_size_position_baseline'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size(latent_dim, window_size),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position'
		), class = type)

	}else if(type == 'vae_fragment_size_position_gmm'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size(latent_dim, window_size),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position'
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_encoder'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position'
		), class = type)

	}else
		stop(sprintf('unknown model type: %s', type))

} # build_model


#' cvae
#'
build_model_cvae <- function(input_dim, feature_dim, hidden_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_cvae(latent_dim, window_size),
		decoder = decoder_model_cvae(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		hidden_dim = hidden_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('cvae'))

} # cvae


#' vae
#'
build_model_vae <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae(latent_dim, window_size),
		decoder = decoder_model_vae(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae'))

} # vae


#' vae_20191216a
#' 
#' This model uses the fragment size vector and poitioin vector as the inputs
#' use the whole v-plot as the output
#'
build_model_vae_20191216a <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae_20191216a(latent_dim, window_size),
		decoder = decoder_model_vae_baseline(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae_20191216a'))

} # vae


#' build_model_vae_20191216b
#' 
#' This model uses the vplot, fragment size vector and poitioin vector as the inputs
#' uses the vplot, fragment size vector and poitioin vector as the outputs as well
#'
build_model_vae_20191216b <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae_20191216b(latent_dim, window_size),
		decoder = decoder_model_vae_20191216b(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae_20191216b'))

} # build_model_vae_20191216b


#' build_model_vae_20191216c
#' 
#' This model uses the vplot, fragment size vector and poitioin vector as the inputs
#' uses the vplot, fragment size vector and poitioin vector as the outputs as well
#'
build_model_vae_20191216c <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae_baseline(latent_dim, window_size),
		decoder = decoder_model_vae_vplot_position_fragment_size(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae_20191216c'))

} # build_model_vae_20191216c


#' build_model_vae_20191216d
#' 
#' This model uses the vplot as the input
#' uses the vplot, and poitioin vector as the outputs as well
#'
build_model_vae_20191216d <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae_baseline(latent_dim, window_size),
		decoder = decoder_model_vae_vplot_position(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae_20191216d'))

} # build_model_vae_20191216d


#' build_model_vae_20191216e
#' 
#' This model uses the vplot as the input and predict and fragment size and poitioin vector 
#'
build_model_vae_20191216e <- function(input_dim, feature_dim, latent_dim, num_samples, window_size){

	structure(list(
		encoder = encoder_model_vae_baseline(latent_dim, window_size),
		decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
		prior = prior_model(latent_dim),
		input_dim = input_dim,
		feature_dim = feature_dim,
		latent_dim = latent_dim,
		num_samples = num_samples,
		window_size = window_size
	), class = c('vae_20191216e'))

} # build_model_vae_20191216e



