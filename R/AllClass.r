setOldClass('kerastools.model.RModel')

setOldClass('tensorflow.tensor')

setClassUnion('listOrNULL', members = c('list', 'NULL'))

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
	'vplot_predictor_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		n_samples = 'integer',
		window_dim = 'integer',
		interval_dim = 'integer',
		latent_dim = 'integer'
	)
)
