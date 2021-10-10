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
#' @import keras 
#' @importFrom tensorflow shape tf_function
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges rowData colData assays rbind rowData<- assays<-
#' @importFrom GenomicRanges resize reduce granges GRangesList slidingWindows width GRanges seqnames start coverage findOverlaps strand
#' @importFrom GenomeInfoDb seqlengths seqlevels seqinfo seqlevels<- seqlengths<- seqinfo<- genome genome<-
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
