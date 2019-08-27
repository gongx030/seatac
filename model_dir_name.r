#' model_dir_name
model_dir_name <- function(dataset, latent_dim, window_size, bin_size, fragment_size_range, fragment_size_interval, epochs){
	f <- sprintf('analysis/seatac/models/dataset=%s_latent_dim=%d_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d_epochs=%d', dataset, latent_dim, window_size, bin_size, fragment_size_range[1], fragment_size_range[2], fragment_size_interval, epochs)
	flog.info(sprintf('model dir: %s', f))
	f
}
