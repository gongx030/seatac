#' model_dir_name
model_dir_name <- function(gs, latent_dim, window_size, bin_size, fragment_size_range, fragment_size_interval, sequence_dim, type, n_components = 0){
	f <- sprintf('analysis/seatac/models/gs=%s_latent_dim=%d_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d_sequence_dim=%d_type=%s_n_components=%d', gs, latent_dim, window_size, bin_size, fragment_size_range[1], fragment_size_range[2], fragment_size_interval, sequence_dim, type, n_components)
	flog.info(sprintf('model dir: %s', f))
	f
}
