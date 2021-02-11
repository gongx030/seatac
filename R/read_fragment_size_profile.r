#' Read fragment size profile
#' 
#' Read fragment size profile from BAM files
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filename BAM file name
#' @param genome a BSgenome object
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(0L, 320L))
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'read_fragment_size_profile',
	signature(
		x = 'GRanges',
		filename = 'character',
		genome = 'BSgenome'
	),
	function(
		x,
		filename, 
		genome,
		fragment_size_range = c(0, 320L),
		fragment_size_interval = 10L
	){	

		stopifnot(length(filename) == 1)

		breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

		g <- read_bam(filename, peaks = x, genome = genome)
		g <- g[strand(g) == '+']
		y <- mcols(g)$isize
		y <- as.numeric(cut(y, breaks))	# discretize the fragment size
		y <- y[!is.na(y)] # remove the read pairs where the fragment size is outside of "fragment_size_range"
		z <- factor(y,  1:length(centers)) %>% table() %>% c()
		class(z) <- 'numeric'
		z <- scale01(z)
		data.frame(
			fragment_size = centers,
			prob = z
		)
	}
)
