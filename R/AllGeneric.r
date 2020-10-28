setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))
setGeneric('read_vplot', function(x, filename, genome, ...) standardGeneric('read_vplot'))
setGeneric('summarize_kmers', function(x, ...) standardGeneric('summarize_kmers'))
setGeneric('count_reads', function(x, filename, genome, ...) standardGeneric('count_reads'))
setGeneric('vplot', function(x, ...) standardGeneric('vplot'))
setGeneric('add_track', function(x, object, ...) standardGeneric('add_track'))
setGeneric('add_kmers', function(x, k, ...) standardGeneric('add_kmers'))
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))
setGeneric('Seq2VplotModel', function(model, ...) standardGeneric('Seq2VplotModel'))
setGeneric('encode', function(model, x, ...) standardGeneric('encode'))
setGeneric('add_kmers_counts', function(x, y, ...) standardGeneric('add_kmers_counts'))
setGeneric('compute_deviations', function(x, annotation, ...) standardGeneric('compute_deviations'))
