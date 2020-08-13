#' seatac: A package for ATAC-seq V-plot analysis
#'
#' ATAC-seq V-plot analysis
#'
#' @import methods
#' @import BiocGenerics
#' @import tfprobability
#' @import tensorflow
#' @import keras
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
#' @importFrom BiocParallel bplapply
#' @import rtracklayer
#' @docType package
#' @name seatac
#'
NULL
# > NULL

