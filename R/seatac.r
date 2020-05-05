#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import tensorflow
#' @import keras
#' @import tfprobability
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians rowMins rowMaxs colVars
#' @import futile.logger 
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments GAlignmentsList
#' @importFrom gplots colorpanel
#' @importFrom abind abind
#' @importFrom fields image.smooth
#' @importFrom purrr pluck
#' @importFrom BSgenome BSgenome
#' @importFrom motifmatchr matchMotifs
#' @importFrom cccd nng
#' @importFrom igraph cluster_louvain membership
#' @docType package
#' @name seatac
#'
NULL
# > NULL

