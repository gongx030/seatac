#' Read BAM file
#'
#' Read the aligned paired-end reads over the given peaks. 
#'
#' @param filename BAM file names 
#' @param peaks a GRange object that define a set of genomic regions.
#' @param genome a BSgenome for genome information
#'
#' @export 
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'  
setMethod(
	'read_bam',
	signature(
		filename = 'character',
		peaks = 'GRanges',
		genome = 'BSgenome'
	), 
	function(filename, peaks, genome){

		validate_bam(filename)

	  # genomic ranges covered by the BAM files
    x <- idxstatsBam(filename)
		gr <- GRanges(seqnames = x[, 'seqnames'], ranges = IRanges(1, x[, 'seqlength']))
		seqlengths(seqinfo(gr)) <- width(gr)
	  genome(seqinfo(gr)) <- metadata(genome)$genome

	  flag <- scanBamFlag(
 	  	isSecondaryAlignment = FALSE,
	    isUnmappedQuery = FALSE,
 	  	isNotPassingQualityControls = FALSE,
	    isProperPair = TRUE
	  )

		# set the seqlevels and seqlengths from BAM files to "which"
		# this block will cause error on R4.1.1. with GenomicAlignments_1.28.0: 
		# Error in .io_bam(.scan_bamfile, file, reverseComplement, yieldSize(file), : seqlevels(param) not in BAM header:

#		seqlevels(peaks, pruning.mode = 'coarse') <- seqlevels(gr)
#	  seqlengths(seqinfo(peaks)) <-  seqlengths(seqinfo(gr))
#	  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))

		# find the seqlevels that are overlapped between BAM and BSgenome
		slevels <- intersect(seqlevels(gr), seqlevels(peaks))
		sinfo <- seqinfo(gr)[slevels]

		# change the seqlevel and seqinfo to the overlapped seqs
		# so that readGAlignments will not produce error when there is missing seqlevels in the BAM
		seqlevels(peaks, pruning.mode = 'coarse') <- slevels
		seqlengths(seqinfo(peaks)) <-  seqlengths(sinfo)
	  genome(seqinfo(peaks)) <-  genome(seqinfo(gr))

 		param <- ScanBamParam(which = reduce(peaks), flag = flag, what = c('isize'), tag = 'RG')

		# Read the PE reads 
		message(sprintf('read_bam | reading %s', filename))
		x <- readGAlignments(filename, param = param)

 		seqlevels(x, pruning.mode = 'coarse') <- seqlevels(gr)
	  seqlengths(seqinfo(x)) <-  seqlengths(seqinfo(gr))
	  genome(seqinfo(x)) <-  genome(seqinfo(gr))

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
#'
#' @author Wuming Gong (gongx030@umn.edu)
#'  
setMethod(
	'read_bam',
	signature(
		filename = 'character',
		peaks = 'missing',
		genome = 'BSgenome'
	), 
	function(filename, peaks, genome){

		validate_bam(filename)

	  # genomic ranges covered by the BAM files
    x <- idxstatsBam(filename)
		gr <- GRanges(seqnames = x[, 'seqnames'], ranges = IRanges(1, x[, 'seqlength']))
		seqlengths(seqinfo(gr)) <- width(gr)
	  genome(seqinfo(gr)) <- metadata(genome)$genome

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

 		seqlevels(x, pruning.mode = 'coarse') <- seqlevels(gr)
	  seqlengths(seqinfo(x)) <-  seqlengths(seqinfo(gr))
	  genome(seqinfo(x)) <-  genome(seqinfo(gr))

		x
	}
) # read_bam

