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
	'vplot_model',
	slot = c(
		window_dim = 'integer',
		interval_dim = 'integer',
		n_samples = 'integer',
		fragment_size_range = 'integer',
		fragment_size_interval = 'integer',
		window_size = 'integer',
		bin_size = 'integer'
	)
)

setClass(
	'vplot_autoencoder_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		latent_dim = 'integer'
	),
	contains = 'vplot_model'
)

setClass(
	'vplot_autoencoder_cluster_model',
	slot = c(
		num_clusters = 'integer',
		sigma = 'numeric',
		gamma = 'numeric',
		membership = 'matrix',
		centers = 'matrix'
	),
	contains = 'vplot_autoencoder_model'
)

setClass(
	'vplot_autoencoder_cluster_v2_model',
	slot = c(
		num_clusters = 'integer',
		sigma = 'numeric',
		gamma = 'numeric',
		membership = 'matrix',
		centers = 'matrix'
	),
	contains = 'vplot_autoencoder_cluster_model'
)

setClass(
	'vplot_autoencoder_rge_model',
	slot = c(
		lambda = 'numeric'
	),
	contains = 'vplot_autoencoder_cluster_model'
)


setClass(
	'vplot_autoencoder_disc_model',
	slot = c(
		lambda = 'numeric'
	),
	contains = 'vplot_autoencoder_cluster_model'
)


