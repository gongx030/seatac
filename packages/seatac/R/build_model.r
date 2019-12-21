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
			output = 'fragment_size_position',
			imputation = FALSE
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
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_encoder_with_imputation'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			imputer = imputer_model_position_fragment_size(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = TRUE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_encoder_gmm'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim, n_components = 10),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_gmm'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size_cnn(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim, n_components = 10),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size_cnn(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_with_imputation'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size_cnn(input_dim, feature_dim, reinterpreted_batch_ndims = 0L),
			prior = prior_model(latent_dim),
			imputer = imputer_model_baseline(input_dim, feature_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = TRUE 
		), class = type)

	}else if(type == 'vae_fragment_size_parametric'){

		 structure(list(
			encoder = encoder_model_vae_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_vplot_parametric(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size',
			output = 'fragment_size',
			imputation = FALSE

		), class = type)

	}else if(type == 'vae_fragment_size_position_mixture_cnn_gmm'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size_mixture_cnn(input_dim, feature_dim),
			prior = prior_model_gmm(latent_dim, n_components = 10),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_mixture_cnn'){

		 structure(list(
			encoder = encoder_model_vae_position_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_position_fragment_size_mixture_cnn(input_dim, feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size'){

		 structure(list(
			encoder = encoder_model_vae_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_fragment_size_cnn(feature_dim),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size',
			output = 'fragment_size',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_gmm'){

		 structure(list(
			encoder = encoder_model_vae_fragment_size_cnn(latent_dim),
			decoder = decoder_model_vae_fragment_size_cnn(feature_dim),
			prior = prior_model_gmm(latent_dim, n_components = 10),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size',
			output = 'fragment_size',
			imputation = FALSE
		), class = type)

	}else if(type == 'vae_fragment_size_position_cnn_kmeans'){

		 structure(list(
			type = type,
			encoder = encoder_model_vae_position_fragment_size_cnn_kmeans(latent_dim, n_components = 5L),
			decoder = decoder_model_vae_position_fragment_size_cnn(input_dim, feature_dim),
			regularizer = regularizer_model_kmeans(latent_dim, n_components = 5L, name = 'kmeans'),
			prior = prior_model(latent_dim),
			input_dim = input_dim,
			feature_dim = feature_dim,
			latent_dim = latent_dim,
			num_samples = num_samples,
			window_size = window_size,
			input = 'fragment_size_position',
			output = 'fragment_size_position',
			imputation = FALSE
		), class = type)

	}else
		stop(sprintf('unknown model type: %s', type))

} # build_model

