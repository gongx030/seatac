#' read_bam
#'
#' The generic function of read_bam
#'
#' @param filename BAM file name
#' @param peaks a peak object
#' @param genome a genome object
#' @param ... Other arguments
#'
setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))

#' read_vplot
#'
#' The generic function of read_vplot
#'
#' @param x a genomic object
#' @param filenames BAM file name
#' @param genome a genome object
#' @param ... Other arguments
#'
setGeneric('read_vplot', function(x, filenames, genome, ...) standardGeneric('read_vplot'))

#' read_aggregated_vplot
#'
#' The generic function of read_aggregated_vplot
#'
#' @param x a genomic object
#' @param filenames BAM file name
#' @param genome a genome object
#' @param ... Other arguments
#'
setGeneric('read_aggregated_vplot', function(x, filenames, genome, ...) standardGeneric('read_aggregated_vplot'))

#' count_reads
#'
#' The generic function of count_reads
#'
#' @param x a genomic object
#' @param filenames BAM file names
#' @param genome a genome object
#' @param ... Other arguments
#'
setGeneric('count_reads', function(x, filenames, genome, ...) standardGeneric('count_reads'))

#' vplot
#'
#' The generic function of vplot
#'
#' @param x a Vplot object
#' @param ... Other arguments
#'
setGeneric('vplot', function(x, ...) standardGeneric('vplot'))

#' prepare_data
#'
#' The generic function of prepare_data
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))

#' fit
#'
#' The generic function of fit
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

#' predict
#'
#' The generic function of predict
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))

#' scale01
#'
#' The generic function of scale01
#'
#' @param x a data object
#' @param ... Other arguments
#'
setGeneric('scale01', function(x, ...) standardGeneric('scale01'))

#' test_accessibility
#'
#' The generic function of test_accessibility
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('test_accessibility', function(model, x, ...) standardGeneric('test_accessibility'))

#' predict_vplots
#'
#' The generic function of predict_vplots
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('predict_vplots', function(model, x, ...) standardGeneric('predict_vplots'))


#' predict_nucleosome
#'
#' The generic function of predict_nucleosome
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('predict_nucleosome', function(model, x, ...) standardGeneric('predict_nucleosome'))


#' read_fragment_size_profile
#'
#' The generic function of read_fragment_size_profile
#'
#' @param x a genomic object
#' @param filename BAM file name
#' @param genome a genome object
#' @param ... Other arguments
#'
setGeneric('read_fragment_size_profile', function(x, filename, genome, ...) standardGeneric('read_fragment_size_profile'))

#' read_fragment_size_profile
#'
#' The generic function of read_fragment_size_profile
#'
#' @param model a Model
#' @param filename Filename of the pretrained model
#' @param ... Other arguments
#'
setGeneric('load_model', function(model, filename, ...) standardGeneric('load_model'))

#' predict_fragment_size
#'
#' The generic function of predict_fragment_size
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('predict_fragment_size', function(model, x, ...) standardGeneric('predict_fragment_size'))

#' add_genomic_features
#'
#' The generic function of add_genomic_features
#'
#' @param x a data object that include the Vplot information
#' @param y a data object 
#' @param ... Other arguments
#'
setGeneric('add_genomic_features', function(x, y, ...) standardGeneric('add_genomic_features'))

#' results
#'
#' The generic function of results
#'
#' @param model a model object
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('results', function(model, x, ...) standardGeneric('results'))

#' call_nucleosomes
#'
#' The generic function of call_nucleosomes
#'
#' @param model a Model
#' @param x a data object that include the Vplot information
#' @param ... Other arguments
#'
setGeneric('call_nucleosomes', function(model, x, ...) standardGeneric('call_nucleosomes'))


#' validate
#'
#' The generic function of validate
#'
#' @param model a Model
#' @param x a data object that include the Vplot information
#'
setGeneric('validate', function(model, x) standardGeneric('validate'))

