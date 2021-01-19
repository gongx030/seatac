setGeneric('read_bam', function(filename, peaks, genome, ...) standardGeneric('read_bam'))
setGeneric('read_vplot', function(x, filenames, genome, ...) standardGeneric('read_vplot'))
setGeneric('count_reads', function(x, filename, genome, ...) standardGeneric('count_reads'))
setGeneric('vplot', function(x, ...) standardGeneric('vplot'))
setGeneric('add_kmers', function(x, k, ...) standardGeneric('add_kmers'))
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))
setGeneric('summarize_vplot', function(x, annotation, ...) standardGeneric('summarize_vplot'))
setGeneric('scale01', function(x, ...) standardGeneric('scale01'))
setGeneric('test_accessibility', function(model, x, ...) standardGeneric('test_accessibility'))
