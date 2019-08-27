
save_windows <- function(windows, gs, window_size, bin_size, fragment_size_range, fragment_size_interval){
	f <- sprintf('%s/windows/dataset=%s_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d.rds', project_dir('seatac'), gs, window_size, bin_size, fragment_size_range[1], fragment_size_range[2], fragment_size_interval)
	flog.info(sprintf('writing %s', f))
	saveRDS(windows, f)
}


load_windows <- function(gs, window_size, bin_size, fragment_size_range, fragment_size_interval){
	f <- sprintf('%s/windows/dataset=%s_window_size=%d_bin_size=%d_fragment_size_range=%d+%d_fragment_size_interval=%d.rds', project_dir('seatac'), gs, window_size, bin_size, fragment_size_range[1], fragment_size_range[2], fragment_size_interval)
	if (!file.exists(f))
		stop(sprintf('%s does not exist', f))
	flog.info(sprintf('reading %s', f))
	windows <- readRDS(f)
	windows
} # load_windows
