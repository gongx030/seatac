setOldClass('kerastools.model.RModel')

setOldClass('tensorflow.tensor')

setOldClass('keras.engine.training.Model')

setClassUnion('listOrNULL', members = c('list', 'NULL'))

setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

setClass(
	'vplot_model',
	slot = c(
		fragment_size_range = 'integer',
		fragment_size_interval = 'integer',
		window_size = 'integer',
		bin_size = 'integer',
		n_bins_per_window = 'integer',
		n_intervals = 'integer',
		n_bins_per_block = 'integer',
		n_blocks_per_window = 'integer',
		block_size = 'integer',
		min_reads_per_block = 'numeric',
		max_reads_per_pixel = 'numeric'
	)
)

setClassUnion('vplot_modelOrNULL', members = c('vplot_model', 'NULL'))

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

setClass(
	'VplotsList',
	contains = 'GRangesList'
)

setClass(
	'sparse_array',
	slot = c(
		subs = 'matrix',
		vals = 'numeric',
		dims = 'numeric',
		dimnames = 'listOrNULL'
	)
)

setClass(
	'sparse_vector',
	slot = c(
		subs = 'numeric',
		vals = 'numeric'
	)
)

setClass(
	'vplot_vae_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		latent_dim = 'integer',
		sigma0 = 'numeric'
	),
	contains = 'vplot_model'
)

