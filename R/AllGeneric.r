#setGeneric('array_reshape', function(x, dims, ...) standardGeneric('array_reshape'))

#setGeneric('get_motif_vplot', function(x, pwms, ...) standardGeneric('get_motif_vplot'))

#setGeneric('as_sparse_array', function(x, ...) standardGeneric('as_sparse_array'))

setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))

setGeneric('get_counts', function(x, filenames, genome, ...) standardGeneric('get_counts'))

#setGeneric('array_permute', function(x, perm, ...) standardGeneric('array_permute'))

#setGeneric('compute_deviations', function(x, ...) standardGeneric('compute_deviations'))

setGeneric('read_vplot', function(x, filenames, genome, ...) standardGeneric('read_vplot'))

setGeneric('add_gc_content', function(x, genome, ...) standardGeneric('add_gc_content'))

setGeneric('build_model', function(name, x, ...) standardGeneric('build_model'))

setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

setGeneric('encode', function(model, x, ...) standardGeneric('encode'))

#setGeneric('as_vplot_autoencoder_cluster_v2_model', function(x, ...) standardGeneric('as_vplot_autoencoder_cluster_v2_model'))

#setGeneric('save_model', function(x, dir, ...) standardGeneric('save_model'))

#setGeneric('load_model', function(dir, ...) standardGeneric('load_model'))

#setGeneric('smooth_vplot', function(x, ...) standardGeneric('smooth_vplot'))

setGeneric('predict', function(model, x, ...) standardGeneric('predict'))

setGeneric('vplot', function(x, ...) standardGeneric('vplot'))

setGeneric('simulate_vplot', function(x, ...) standardGeneric('simulate_vplot'))

setGeneric('prepare_vplot', function(x, ...) standardGeneric('prepare_vplot'))

#setGeneric('results', function(x, targets, ...) standardGeneric('results'))

setGeneric('extract_blocks', function(x, ...) standardGeneric('extract_blocks'))

setGeneric('add_nucleosome_signal', function(x, ...) standardGeneric('add_nucleosome_signal'))

setGeneric('add_track', function(x, ...) standardGeneric('add_track'))

setGeneric('add_kmers', function(x, k, genome, ...) standardGeneric('add_kmers'))
