setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))
setGeneric('read_vplot', function(x, filenames, genome, ...) standardGeneric('read_vplot'))
setGeneric('build_model', function(name, x, ...) standardGeneric('build_model'))

setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

setGeneric('encode', function(model, x, ...) standardGeneric('encode'))

#setGeneric('as_vplot_autoencoder_cluster_v2_model', function(x, ...) standardGeneric('as_vplot_autoencoder_cluster_v2_model'))

#setGeneric('save_model', function(x, dir, ...) standardGeneric('save_model'))

#setGeneric('load_model', function(dir, ...) standardGeneric('load_model'))

#setGeneric('smooth_vplot', function(x, ...) standardGeneric('smooth_vplot'))

setGeneric('predict', function(model, x, ...) standardGeneric('predict'))

setGeneric('evaluate', function(model, x, ...) standardGeneric('evaluate'))

setGeneric('predict_nucleosome', function(model, x, ...) standardGeneric('predict_nucleosome'))

setGeneric('vplot', function(x, ...) standardGeneric('vplot'))

#setGeneric('results', function(x, targets, ...) standardGeneric('results'))

#setGeneric('extract_blocks', function(x, ...) standardGeneric('extract_blocks'))

setGeneric('add_track', function(x, ...) standardGeneric('add_track'))

setGeneric('add_kmers', function(x, k, genome, ...) standardGeneric('add_kmers'))

setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))

