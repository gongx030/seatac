setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))
setGeneric('read_vplot', function(x, filename, genome, ...) standardGeneric('read_vplot'))
setGeneric('count_reads', function(x, filename, genome, ...) standardGeneric('count_reads'))
setGeneric('vplot', function(x, ...) standardGeneric('vplot'))
setGeneric('add_track', function(x, object, ...) standardGeneric('add_track'))
setGeneric('add_kmers', function(x, k, ...) standardGeneric('add_kmers'))
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))
setGeneric('Seq2VplotModel', function(model, ...) standardGeneric('Seq2VplotModel'))
setGeneric('encode', function(model, x, ...) standardGeneric('encode'))
setGeneric('decode', function(model, x, ...) standardGeneric('decode'))
setGeneric('summarize_vplot', function(x, annotation, ...) standardGeneric('summarize_vplot'))
setGeneric('scale01', function(x, ...) standardGeneric('scale01'))
setGeneric('cluster_vplot', function(x, ...) standardGeneric('cluster_vplot'))
setGeneric('select_vplot', function(x, ...) standardGeneric('select_vplot'))
setGeneric('match_motifs', function(x, pwms, ...) standardGeneric('match_motifs'))
