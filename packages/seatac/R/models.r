
#' save_model
#'
save_model <- function(x, dir){

	if (missing(dir))
		stop('dir must be specified')

	if (!file.exists(dir))
		dir.create(dir, recursive = TRUE)

	vae_file <- sprintf('%s/vae.h5', dir)
	flog.info(sprintf('writing %s', vae_file))
	x$vae$save_weights(vae_file)
	x$vae_file <- vae_file
	x$vae <- NULL

	latent_file <- sprintf('%s/latent.h5', dir)
	flog.info(sprintf('writing %s', latent_file))
	x$latent$save_weights(latent_file)
	x$latent_file <- latent_file
	x$latent <- NULL

	sequence_motifs_file <- sprintf('%s/sequence_motifs.h5', dir)
	flog.info(sprintf('writing %s', sequence_motifs_file))
	x$sequence_motifs$save_weights(sequence_motifs_file)
	x$sequence_motifs_file <- sequence_motifs_file
	x$sequence_motifs <- NULL

	model_file <- sprintf('%s/model.rds', dir)
	flog.info(sprintf('writing %s', model_file))
	saveRDS(x, model_file)

} # save_model


#' load_model
load_model <- function(dir){

	if (missing(dir))
		stop('dir must be specified')

	model_file <- sprintf('%s/model.rds', dir)

	if (!file.exists(model_file))
		stop(sprintf('%s does not exist', model_file))

	flog.info(sprintf('reading %s', model_file))
	x <- readRDS(model_file)

	model <- cvae(
		input_dim = x$input_dim,
		feature_dim = x$feature_dim,
		latent_dim = x$latent_dim,
		num_samples = x$num_samples,
		window_size = x$window_size ,
		sequence_dim = x$sequence_dim,
		is_nfr = x$is_nfr,
		is_mono_nucleosome = x$is_mono_nucleosome
	)

	flog.info(sprintf('reading %s', x$vae_file))
	model$vae$load_weights(x$vae_file)

	flog.info(sprintf('reading %s', x$latent_file))
	model$latent$load_weights(x$latent_file)

	flog.info(sprintf('reading %s', x$sequence_motifs_file))
	model$sequence_motifs$load_weights(x$sequence_motifs_file)

	model

} # load_model

