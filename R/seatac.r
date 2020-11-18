#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import reticulate
#' @import methods
#' @import Matrix
#' @import dplyr
#' @import tfdatasets
#' @importFrom keras keras_model_custom save_model_weights_tf load_model_weights_tf
#' @importFrom tensorflow shape
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays
#' @importFrom GenomicRanges resize reduce granges
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments GAlignmentsList
#' @importFrom gplots colorpanel
#' @importFrom BSgenome BSgenome
#' @importFrom RANN nn2
#' @importFrom igraph graph.adjacency cluster_louvain
#' @importFrom abind abind
#' @docType package
#' @name seatac
#'
NULL
# > NULL
