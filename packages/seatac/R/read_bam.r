#' read_bam
#'
read_bam <- function(
	filenames, 
	peaks, 
	genome, 
	expand = 2000
){

	validate_bam(filenames)

  num_samples <- length(filenames)

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
  }))
  seqlengths(seqinfo(gr)) <- width(gr)
  genome(seqinfo(gr)) <- providerVersion(genome)

	# set the seqlevels and seqlengths from BAM files to "which"
  seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(gr)
  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(peaks)
  seqlengths(seqinfo(peaks)) <-  seqlengths(seqinfo(gr))
  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))

	flog.info(sprintf('expanding each peak to [-%d, +%d] region', expand / 2, expand / 2))
	windows <- resize(peaks, fix = 'center', width = expand)	# expand the input peaks

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )
  param <- ScanBamParam(which = reduce(windows), flag = flag, what = 'isize')

	# Read the PE reads overlapping with specified windows
  x <- Reduce('c', lapply(1:num_samples, function(i){
    flog.info(sprintf('reading %s', filenames[i]))
    ga <- readGAlignments(filenames[i], param = param)
		mcols(ga)$group <- i
		ga
	}))
	metadata(x)$num_samples <- num_samples
	x

} # readBAM


