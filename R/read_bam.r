#' read_bam
#'
#' Read the aligned paired-end reads over the given peaks. 
#'  
read_bam <- function(
	filenames, 
	peaks, 
	genome
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

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isProperPair = TRUE
  )

	if (!missing(peaks)){
		# set the seqlevels and seqlengths from BAM files to "which"
 		seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(gr)
	  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(peaks)
	  seqlengths(seqinfo(peaks)) <-  seqlengths(seqinfo(gr))
	  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))
  	param <- ScanBamParam(which = reduce(peaks), flag = flag, what = 'isize')
	}else{
  	param <- ScanBamParam(flag = flag, what = 'isize')
	}

	# Read the PE reads 
  x <- Reduce('c', lapply(1:num_samples, function(i){
    flog.info(sprintf('reading %s', filenames[i]))
    ga <- readGAlignments(filenames[i], param = param)
		ga
	}))
	metadata(x)$num_samples <- num_samples
	x

} # read_bam


