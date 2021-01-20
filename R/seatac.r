#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import reticulate
#' @import methods
#' @import Matrix
#' @import dplyr
#' @import tfdatasets
#' @import tfprobability
#' @importFrom keras keras_model_custom save_model_weights_tf load_model_weights_tf
#' @importFrom tensorflow shape
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays rbind
#' @importFrom GenomicRanges resize reduce granges GRangesList slidingWindows
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments GAlignmentsList
#' @importFrom gplots colorpanel
#' @importFrom BSgenome BSgenome
#' @docType package
#' @name seatac
#'
NULL
# > NULL
