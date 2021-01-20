#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filenames BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(80L, 320L))
#' @param fragment_size_interval fragment_size_interval
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'read_vplot',
	signature(
		x = 'GRanges',
		filenames = 'character',
		genome = 'BSgenome'
	), 
	function(
		x, 
		filenames, 
		genome,
		bin_size = 5L,
		fragment_size_range = c(80, 320),
		fragment_size_interval = 5L
	){

		if (is.null(names(filenames)))
			names(filenames) <- 1:length(filenames)

		se <- lapply(1:length(filenames), function(i){
			xi <- read_vplot_core(x, filenames[i], genome, bin_size, fragment_size_range, fragment_size_interval)
			rowData(xi)$batch <- names(filenames)[i]
			xi@samples <- names(filenames)[i]
			xi
		})
		se <- do.call('rbind', se)
		se@samples <- unique(rowData(se)$batch)
		se@n_samples <- length(se@samples)
		se
	}

) # read_vplot


#' Read aggregated V-plot from a list of genomic ranges (GRangeList)
#'
#' This function is to quickly read aggregated V-plot given a list of GRanges (e.g. motif binding sites) for 
#' scATAC-seq.
#' 
#' @param x a GRangeList object that define a list of genomic region sets.
#' @param filenames BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size in base pairs (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes (default: c(80L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 5L)
#' @return a SummarizedVplots object
#'
#' @export
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'read_vplot',
	signature(
		x = 'GRangesList',
		filenames = 'character',
		genome = 'BSgenome'
	), 
	function(
		x, 
		filenames, 
		genome,
		bin_size = 5L,
		fragment_size_range = c(80, 320),
		fragment_size_interval = 5L
	){

		stopifnot(length(filenames) == 1)

		window_size <- width(x) %>% unlist()

		if (length(unique(window_size)) > 1)
			stop('the window size of input data must be equal')

		window_size <- window_size[1]

		stopifnot(window_size %% bin_size == 0)

 		n_bins_per_window <- window_size / bin_size

		breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

		n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval

		classes <- names(x)	# could be no names
		if (is.null(classes))
			classes <- 1:length(x)
		idx <- rep(1:length(x), sapply(x, length))	# number of binding site for each group
		x <- unlist(x)
		names(x) <- NULL

		# read everything from the BAM file.  For bulk ATAC-seq, this may change to only read a subset of regions
		g <- read_bam(filenames, genome = genome)
		g <- g[strand(g) == '+']
		g <- GRanges(
			seqnames = seqnames(g), 
			range = IRanges(start(g) + round(mcols(g)$isize / 2), width = 1), 
			isize = mcols(g)$isize
		)
		g$fragment_size <- as.numeric(cut(g$isize, breaks))	# discretize the fragment size
		g <- g[!is.na(g$fragment_size)] # remove the read pairs where the fragment size is outside of "fragment_size_range"

		if (length(g) == 0){
			warning(sprintf('there is no qualified reads in %s', filenames))
			return(NULL)
		}

		non_empty <- x %over% g	

		if (!any(non_empty)){
			sprintf('there is no qualified motifs in %s', filenames) %>% warning()
			return(NULL)
		}

		x <- x[non_empty]
		idx <- idx[non_empty]

		MW <- sparseMatrix(i = idx, j = 1:length(x), dims = c(length(classes), length(x)), dimnames = list(classes, NULL))	# motif ~ windows

		bins <- x %>%
			slidingWindows(width = bin_size, step = bin_size) %>%
			unlist()

		CF <- sparseMatrix(i = 1:length(g), j = g$fragment_size, dims = c(length(g), n_intervals))  # read center ~ fragment size
		BC <- as.matrix(findOverlaps(bins, g))	# bins ~ read center
		BC <- as(sparseMatrix(BC[, 1], BC[, 2], dims = c(length(bins), length(g))), 'dgCMatrix') # bins ~ read center

		# a transpose operator for converting [n_intervals, n_bins_per_window] to [n_bins_per_window, n_intervals]
		h <- c(t(matrix(1:c(n_bins_per_window * n_intervals), n_intervals, n_bins_per_window)))	

		counts <- t(BC %*% CF)  # bins ~ fragment size
		dim(counts) <- c(n_bins_per_window * n_intervals, length(x))
		counts <- t(counts)
		counts <- counts[, h]

		counts <- MW %*% counts	# motif ~ n_bins_per_window * n_intervals

		se <- SummarizedExperiment(
			assays = list(counts = counts)
		)

		new(
			'SummarizedVplots', 
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
)



#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filename BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(80L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 5L)
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
read_vplot_core <- function(
	x, 
	filename, 
	genome, 
	bin_size = 5L, 
	fragment_size_range = c(80L, 320L), 
	fragment_size_interval = 5L
){
	window_size <- width(x)

	if (length(unique(window_size)) > 1)
		stop('the window size of input data must be equal')

	window_size <- window_size[1]

	stopifnot(window_size %% bin_size == 0)

	# compute the center point between PE reads
	# this is faster than using GAlignmentPairs
	window_size <- width(x)[1]
	n_bins_per_window <- window_size / bin_size
	breaks <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
	centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
	n_intervals <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval
	bins <- x %>%
		slidingWindows(width = bin_size, step = bin_size) %>%
		unlist()
	n_bins_per_window <- window_size / bin_size
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
	rowData(se)$id <- 1:length(x)

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
		positions = seq(bin_size, window_size, by = bin_size) - (window_size / 2),
		n_samples = 1L
	)
} # read_vplot_core

