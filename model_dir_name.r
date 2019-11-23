#' model_dir_name
model_dir_name <- function(gs, latent_dim, window_size, bin_size, fragment_size_range, fragment_size_interval, sequence_dim, n_components = 0, min_reads_per_window = 0, flanking, epochs, batch_size, steps_per_epoch){
	f <- sprintf('analysis/seatac/models/gs=%s_latent_dim=%d_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d_sequence_dim=%d_n_components=%d_min_reads_per_window=%d_flanking=%d_epochs=%d_batch_size=%d_steps_per_epoch=%d', 
		gs, 
		latent_dim, 
		window_size, 
		bin_size, 
		fragment_size_range[1], 
		fragment_size_range[2], 
		fragment_size_interval, 
		sequence_dim, 
		n_components, 
		min_reads_per_window,
		flanking,
		epochs,
		batch_size,
		steps_per_epoch
	)
	flog.info(sprintf('model dir: %s', f))
	f
} # model_dir_name
