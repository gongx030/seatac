#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import reticulate
#' @import methods
#' @import Matrix
#' @import futile.logger 
#' @import dplyr
#' @import rhdf5
#' @importFrom keras keras_model_custom
#' @importFrom tensorflow shape
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays
#' @importFrom GenomicRanges resize reduce granges
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
