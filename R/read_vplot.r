#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filenames BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(0, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param ignore_strand whether ignore the strand of the V-plot (default: TRUE)
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
		fragment_size_range = c(0, 320L),
		fragment_size_interval = 10L,
		ignore_strand = TRUE
	){

		if (is.null(names(filenames)))
			names(filenames) <- 1:length(filenames)

		se <- lapply(1:length(filenames), function(i){
			xi <- read_vplot_core(x, filenames[i], genome, bin_size, fragment_size_range, fragment_size_interval, ignore_strand = ignore_strand)
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


#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRangesList object that define a set of genomic regions.
#' @param filenames BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(0, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param ignore_strand whether ignore the strand of the V-plot (default: TRUE)
#'
#' @export
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
		fragment_size_range = c(0, 320L),
		fragment_size_interval = 10L,
		ignore_strand = TRUE
	){

		if (is.null(names(filenames)))
			names(filenames) <- 1:length(filenames)

		stopifnot(!any(duplicated(names(filenames))))

		stopifnot(length(x) == length(filenames))

		se <- lapply(1:length(filenames), function(i){
			xi <- read_vplot_core(x[[i]], filenames[i], genome, bin_size, fragment_size_range, fragment_size_interval, ignore_strand = ignore_strand)
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



#' Read the V-plot
#' 
#' Read the V-plot from BAM files within a set of genomic regions
#'
#' @param x a GRange object that define a set of genomic regions.
#' @param filename BAM file names
#' @param genome a BSgenome object
#' @param bin_size The bin size (default: 5L)
#' @param fragment_size_range The range of the PE reads fragment sizes that are used for 
#'				constructing Vplot (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param ignore_strand whether ignore the strand of the V-plot (default: TRUE)
#' @importFrom GenomicRanges start end
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'
read_vplot_core <- function(
	x, 
	filename, 
	genome, 
	bin_size = 5L, 
	fragment_size_range = c(0L, 320L), 
	fragment_size_interval = 10L,
	ignore_strand = TRUE
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
	n_bins_per_window <- window_size / bin_size

	starts <- matrix(start(x), nrow = n_bins_per_window, ncol = length(x), byrow = TRUE)
	starts <- starts + cumsum(rep(bin_size, n_bins_per_window)) - bin_size

	if (!ignore_strand){

		if (any(as.character(strand(x)) %in% '*'))
			stop('when ignore_strand=TRUE, strand(x) must be either "+" or "-"')

		is_minus <- as.character(strand(x)) == '-'

		if (any(is_minus)){
			starts[, is_minus] <- matrix(end(x)[is_minus], nrow = n_bins_per_window, ncol = sum(is_minus) , byrow = TRUE)
			starts[, is_minus] <- starts[, is_minus] - cumsum(rep(bin_size, n_bins_per_window)) + 1L
		}
	}

	bins <- GRanges(seqnames = rep(seqnames(x), each = n_bins_per_window), ranges = IRanges(start = c(starts), width = bin_size))

	peaks <- reduce(resize(x, fix = 'center', width = window_size + 2000))
	g <- read_bam(filename, peaks = peaks, genome = genome)

	g <- g[strand(g) == '+']
	g <- GRanges(
		seqnames = seqnames(g), 
		ranges = IRanges(start(g) + round(mcols(g)$isize / 2), width = 1), 
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
	counts <- counts[, h, drop = FALSE]

	seqlevels(x) <- seqlevels(genome)
	seqlengths(seqinfo(x)) <- seqlengths(genome)
	genome(seqinfo(x)) <- metadata(genome)$genome

	se <- SummarizedExperiment(
		assays = list(counts = counts),
		rowRanges = x
	)
	rowData(se)$id <- 1:length(x)
	rowData(se)$sequence <- getSeq(genome, x)

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

