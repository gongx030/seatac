#' read_bam
#'
#' The generic function of read_bam
#'
#' @param filename BAM file name
#' @param peaks a peak object
#' @param genome a genome object
#'
setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))

#' read_vplot
#'
#' The generic function of read_vplot
#'
#' @param filename BAM file name
#' @param peaks a peak object
#' @param genome a genome object
#'
setGeneric('read_vplot', function(x, filenames, genome, ...) standardGeneric('read_vplot'))

#' count_reads
#'
#' The generic function of count_reads
#'
#' @param filename BAM file name
#' @param peaks a peak object
#' @param genome a genome object
#'
setGeneric('count_reads', function(x, filename, genome, ...) standardGeneric('count_reads'))

#' vplot
#'
#' The generic function of vplot
#'
#' @param x a Vplot object
#'
setGeneric('vplot', function(x, ...) standardGeneric('vplot'))

#' prepare_data
#'
#' The generic function of prepare_data
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#'
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))

#' fit
#'
#' The generic function of fit
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#'
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

#' predict
#'
#' The generic function of predict
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#'
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))

#' scale01
#'
#' The generic function of scale01
#'
#' @param x a data object
#'
setGeneric('scale01', function(x, ...) standardGeneric('scale01'))

#' test_accessibility
#'
#' The generic function of test_accessibility
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#'
setGeneric('test_accessibility', function(model, x, ...) standardGeneric('test_accessibility'))
