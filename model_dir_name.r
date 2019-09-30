#' model_dir_name
model_dir_name <- function(gs, latent_dim, window_size, bin_size, fragment_size_range, fragment_size_interval, sequence_dim){
	f <- sprintf('analysis/seatac/models/gs=%s_latent_dim=%d_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d_sequence_dim=%d', gs, latent_dim, window_size, bin_size, fragment_size_range[1], fragment_size_range[2], fragment_size_interval, sequence_dim)
	flog.info(sprintf('model dir: %s', f))
	f
}
