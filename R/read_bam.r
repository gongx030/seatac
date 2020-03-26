setGeneric('read_bam', function(filenames, peaks, genome, ...) standardGeneric('read_bam'))

#' read_bam
#'
#' Read the aligned paired-end reads over the given peaks. 
#' @param filenames BAM file names 
#' @param peaks a GRange object that define a set of genomic regions.
#' @param genome a BSgenome for genome information
#'  
setMethod(
	'read_bam',
	signature(
		filenames = 'character',
		peaks = 'GRanges',
		genome = 'BSgenome'
	), 
	function(filenames, peaks, genome, ...){

		validate_bam(filenames)

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

		# set the seqlevels and seqlengths from BAM files to "which"
 		seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(gr)
	  seqlevels(gr, pruning.mode = 'coarse') <- seqlevels(peaks)
	  seqlengths(seqinfo(peaks)) <-  seqlengths(seqinfo(gr))
	  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))
 		param <- ScanBamParam(which = reduce(peaks), flag = flag, what = 'isize')

		# Read the PE reads 
	  x <- lapply(seq_len(length(filenames)), function(i){
 	  	flog.info(sprintf('reading %s', filenames[i]))
	    readGAlignments(filenames[i], param = param)
		}) %>% 
			GAlignmentsList()
		x
	}
) # read_bam


setMethod(
	'read_bam',
	signature(
		filenames = 'character',
		peaks = 'missing',
		genome = 'BSgenome'
	), 
	function(filenames, peaks, genome, ...){

		read_bam(filenames)

	}
) # read_bam


setMethod(
	'read_bam',
	signature(
		filenames = 'character',
		peaks = 'missing',
		genome = 'missing'
	), 
	function(filenames, peaks, genome, ...){

		validate_bam(filenames)

	  flag <- scanBamFlag(
 	  	isSecondaryAlignment = FALSE,
	    isUnmappedQuery = FALSE,
 	  	isNotPassingQualityControls = FALSE,
	    isProperPair = TRUE
	  )

 		param <- ScanBamParam(flag = flag, what = 'isize')

		# Read the PE reads 
	  x <- lapply(seq_len(length(filenames)), function(i){
 	  	flog.info(sprintf('reading %s', filenames[i]))
	    readGAlignments(filenames[i], param = param)
		}) %>% 
			GAlignmentsList()
		x
	}
) # read_bam

