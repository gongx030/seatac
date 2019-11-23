#' make_windows
#'
make_windows <- function(x, bam_file, window_size = 640, bin_size = 5, fragment_size_range = c(50, 370), fragment_size_interval = 5, genome){

	x <- resize(x, fix = 'center', width = window_size)
	mcols(x)$sequence <- getSeq(genome, x)

	ga <- read_bam(bam_file, x, genome = genome, expand = window_size * 2)
	x <- readFragmentSizeMatrix(ga, x, window_size = window_size, bin_size = bin_size, fragment_size_range = fragment_size_range, fragment_size_interval = fragment_size_interval)
	x

} # make_windows
