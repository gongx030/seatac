setOldClass('tensorflow.tensor')

setOldClass('tf_dataset')

setClassUnion('listOrNULL', members = c('list', 'NULL'))

setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

setClass(
	'Vplots', 
	slot = c(
		fragment_size_range  = 'integer',
		fragment_size_interval = 'integer',
		bin_size = 'integer',
		window_size = 'integer',
		n_intervals = 'integer',
		n_bins_per_window = 'integer',
		breaks = 'numeric',
		centers = 'numeric',
		positions = 'numeric'
	),
	contains = 'GRanges'
)

setClass(
	'VplotsKmers',
	slot = c(
		kmers = 'character',
		k = 'integer'
	),
	contains = 'Vplots'
)

