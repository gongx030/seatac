#' validate_bam
#'
#' Validate BAM files
#'
#' @param filenames BAM file names
#'
validate_bam <- function(filenames){

  if (missing(filenames))
    stop('filenames are missing')

  existed <- file.exists(filenames)
  if (any(!existed))
    stop(sprintf('%s do not exist', paste(filenames[!existed], collapse = ', ')))

	# check if the BAM index file exists
	index_files <- sprintf('%s.bai', filenames)
  existed <- file.exists(index_files)
  if (any(!existed)){
		message(sprintf('validate_bam | indexing bam files: %s', paste(filenames[!existed], collapse = ',')))
		indexBam(filenames[!existed])
	}

  is_pe <- sapply(filenames, testPairedEndBam)
  if(any(!is_pe)){
    stop(paste(filenames[!s_pe], collapse = ', '),"are not paired-end files.")
  }

} # validate_bam
