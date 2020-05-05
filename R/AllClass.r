setOldClass('kerastools.model.RModel')

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
	'vplot_autoencoder_model',
	slot = c(
		batch_correcter = 'kerastools.model.RModel',
		encoder = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		latent_dim = 'integer',
		n_samples = 'integer',
		window_dim = 'integer',
		interval_dim = 'integer'
	)
)
