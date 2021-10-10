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
#' @importFrom tensorflow shape tf_function
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays rbind rowData<- assays<-
#' @import GenomicRanges 
#' @import GenomeInfoDb 
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam indexBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments GAlignmentsList
#' @importFrom gplots colorpanel
#' @importFrom BSgenome BSgenome getSeq
#' @importFrom S4Vectors metadata mcols mcols<-
#' @importFrom abind abind
#' @importFrom stats pchisq p.adjust
#' @importFrom IRanges IRanges %over%
#' @importFrom graphics abline axis
#' @docType package
#' @name seatac
#'
NULL
# > NULL
