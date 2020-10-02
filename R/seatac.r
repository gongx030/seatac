#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import methods
#' @import Matrix
#' @import futile.logger 
#' @import dplyr
#' @importFrom GenomicRanges resize reduce
#' @importFrom reticulate PyClass
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments GAlignmentsList
#' @importFrom gplots colorpanel
#' @importFrom BSgenome BSgenome
#' @importFrom BiocParallel bplapply
#' @docType package
#' @name seatac
#'
NULL
# > NULL
