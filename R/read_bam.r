#' Read BAM file
#'
#' Read the aligned paired-end reads over the given peaks. 
#'
#' @param filename BAM file names 
#' @param peaks a GRange object that define a set of genomic regions.
#' @param genome a BSgenome for genome information
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'  
setMethod(
	'read_bam',
	signature(
		filename = 'character',
		peaks = 'GRanges',
		genome = 'BSgenome'
	), 
	function(filename, peaks, genome, ...){

		validate_bam(filename)

	  # genomic ranges covered by the BAM files
    x <- idxstatsBam(filename)
		gr <- GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
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
		message(sprintf('read_bam | reading %s', filename))
		x <- readGAlignments(filename, param = param)
		x
	}
) # read_bam



#' Read BAM file
#'
#' Read all aligned paired-end reads
#'
#' @param filename BAM file names 
#' @param peaks missing
#' @param genome a BSgenome for genome information
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'  
setMethod(
	'read_bam',
	signature(
		filename = 'character',
		peaks = 'missing',
		genome = 'BSgenome'
	), 
	function(filename, peaks, genome, ...){

		validate_bam(filename)

	  # genomic ranges covered by the BAM files
    x <- idxstatsBam(filename)
		gr <- GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
		seqlengths(seqinfo(gr)) <- width(gr)
	  genome(seqinfo(gr)) <- providerVersion(genome)

	  flag <- scanBamFlag(
 	  	isSecondaryAlignment = FALSE,
	    isUnmappedQuery = FALSE,
 	  	isNotPassingQualityControls = FALSE,
	    isProperPair = TRUE
	  )

 		param <- ScanBamParam(flag = flag, what = 'isize')

		# Read the PE reads 
		message(sprintf('read_bam | reading %s', filename))
		x <- readGAlignments(filename, param = param)
		x
	}
) # read_bam

