#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filename BAM file name
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes (default: c(80L, 320L))
#' @param fragment_size_interval fragment_size_interval
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'read_vplot',
	signature(
		x = 'GRanges',
		filename = 'character',
		genome = 'BSgenome'
	), 
	function(
		x, 
		filename, 
		genome,
		bin_size = 5L,
		fragment_size_range = c(80, 320),
		fragment_size_interval = 5L
	){

		window_size <- width(x)

		if (length(unique(window_size)) > 1)
			stop('the window size of input data must be equal')

		window_size <- window_size[1]

		stopifnot(window_size %% bin_size == 0)

 		n_bins_per_window <- window_size / bin_size

		breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

		n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval

		bins <- x %>%
			slidingWindows(width = bin_size, step = bin_size) %>%
			unlist()

		# compute the center point between PE reads
		# this is faster than using GAlignmentPairs
		peaks <- reduce(resize(x, fix = 'center', width = window_size + 2000))

		g <- read_bam(filename, peaks = peaks, genome = genome)
		g <- g[strand(g) == '+']
		g <- GRanges(
			seqnames = seqnames(g), 
			range = IRanges(start(g) + round(mcols(g)$isize / 2), width = 1), 
			isize = mcols(g)$isize
		)
		g$fragment_size <- as.numeric(cut(g$isize, breaks))	# discretize the fragment size
		g <- g[!is.na(g$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"

		if (length(g) == 0){
			warning(sprintf('there is no qualified reads in %s', filename))
			return(NULL)
		}

		CF <- sparseMatrix(i = 1:length(g), j = g$fragment_size, dims = c(length(g), n_intervals))  # read center ~ fragment size
		BC <- as.matrix(findOverlaps(bins, g))	# bins ~ read center
		BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(g))), 'dgCMatrix') # bins ~ read center

		# a transpose operator for converting [n_intervals, n_bins_per_window] to [n_bins_per_window, n_intervals]
		h <- c(t(matrix(1:c(n_bins_per_window * n_intervals), n_intervals, n_bins_per_window)))	

		counts <- t(BC %*% CF)  # bins ~ fragment size
		dim(counts) <- c(n_bins_per_window * n_intervals, length(x))
		counts <- t(counts)
		counts <- counts[, h]

		seqlevels(x) <- seqlevels(genome)
		seqlengths(seqinfo(x)) <- seqlengths(genome)
		genome(seqinfo(x)) <- providerVersion(genome)

		se <- SummarizedExperiment(
			assays = list(counts = counts),
			rowRanges = x
		)

		new(
			'Vplots', 
			se, 
			fragment_size_range  = as.integer(fragment_size_range),
			fragment_size_interval = as.integer(fragment_size_interval),
			bin_size = as.integer(bin_size),
			window_size = as.integer(window_size),
			n_intervals = as.integer(n_intervals),
			n_bins_per_window = as.integer(n_bins_per_window ),
			breaks = breaks,
			centers = centers,
			positions = seq(bin_size, window_size, by = bin_size) - (window_size / 2)
		)
	}

) # read_vplot


#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filename BAM file name
#' @param genome genome abbreviation such as mm10 or hg19
#' @param ...
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'read_vplot',
	signature(
		x = 'GRanges',
		filename = 'character',
		genome = 'character'
	), 
	function(
		x, 
		filename, 
		genome,
		...
	){
		genome <- get_bsgenome(genome)
		read_vplot(x, filename, genome, ...)
	}
)
