#' Count reads
#'
#' Count how many reads fall into a specific fragment size range
#'
#' @param x a GRanges object defining the peaks; the width must be the same for all peaks.
#' @param filenames BAM file name(s).
#' @param genome a BS genome object such as BSgenome.Mmusculus.UCSC.mm10
#' @param fragment_size_range fragment size ranges (default:  c(80, 320))
#' @return read counts
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'count_reads',
	signature(
		x = 'GRanges',
		filenames = 'character',
		genome = 'BSgenome'
	), 
	function(
		x, 
		filenames, 
		genome,
		fragment_size_range = c(0L, 320)
	){

		if (is.null(names(filenames)))
			names(filenames) <- 1:length(filenames)

		do.call('cbind', lapply(1:length(filenames), function(i){
			count_reads_core(x, filenames[i], genome, fragment_size_range)
		}))

	}
) # count_reads



#' Count reads
#'
#' Count how many reads fall into a specific fragment size range
#'
#' @param x a GRangesList object defining the peaks; the width must be the same for all peaks.
#' @param filenames BAM file name(s).
#' @param genome a BS genome object such as BSgenome.Mmusculus.UCSC.mm10
#' @param fragment_size_range fragment size ranges (default:  c(80, 320))
#' @return a list of read counts
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'count_reads',
	signature(
		x = 'GRangesList',
		filenames = 'character',
		genome = 'BSgenome'
	), 
	function(
		x, 
		filenames, 
		genome,
		fragment_size_range = c(0L, 320)
	){

		if (length(x) != length(filenames))
			stop('The length of x should be the same as filenames')

		counts <- lapply(1:length(filenames), function(i){
			count_reads(x[[i]], filenames[i], genome)
		})
		counts
	}
) # count_reads

	

#' count_reads_core
#'
#' Count how many reads fall into a specific fragment size range
#'
#' @param x a GRanges object defining the peaks; the width must be the same for all peaks.
#' @param filename BAM file name.
#' @param genome a BS genome object such as BSgenome.Mmusculus.UCSC.mm10
#' @param fragment_size_range fragment size ranges (default:  c(80, 320))
#' @importFrom GenomicRanges trim
#'
#' @return read counts
#' 
#' @author Wuming Gong (gongx030@umn.edu)
#'
count_reads_core <- function(
	x, 
	filename, 
	genome, 
	fragment_size_range = c(0L, 320)
){

	window_size <- width(x)

	if (length(unique(window_size)) > 1)
		stop('the window size of input data must be equal')

	peaks <- reduce(resize(x, fix = 'center', width = window_size * 2))

	g <- read_bam(filename, peaks = peaks, genome = genome)

	seqlevels(x, pruning.mode = 'coarse') <- seqlevels(g)
  seqlengths(seqinfo(x)) <-  seqlengths(seqinfo(g))
 	genome(seqinfo(x)) <-  genome(seqinfo(g))
	x <- trim(x)

	g <- g[strand(g) == '+']
	g <- GRanges(
		seqnames = seqnames(g), 
		ranges = IRanges(start(g) + round(mcols(g)$isize / 2), width = 1), 
		isize = mcols(g)$isize,
		seqinfo = seqinfo(g),
		seqlengths = seqlengths(g)
	)
	g <- g[g$isize >= fragment_size_range[1] & g$isize <= fragment_size_range[2]]

	coverage(g)[x] %>% sum()
}

