setOldClass('kerastools.model.RModel')

setOldClass('tensorflow.tensor')

setClassUnion('listOrNULL', members = c('list', 'NULL'))

setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))

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
		n_samples = 'integer',
		fragment_size_range = 'integer',
		fragment_size_interval = 'integer',
		window_size = 'integer',
		bin_size = 'integer',
		n_bins_per_window = 'integer',
		n_intervals = 'integer'
	)
)

setClass(
	'vplot_parametric_vae_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder_window = 'kerastools.model.RModel',
		decoder_interval = 'kerastools.model.RModel',
		latent_dim = 'integer',
		sigma0 = 'numeric'
	),
	contains = 'vplot_model'
)

setClass(
	'vplot_parametric_vae_v2_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		latent_dim = 'integer',
		sigma0 = 'numeric'
	),
	contains = 'vplot_model'
)

setClass(
	'vplot_parametric_vae_v3_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		latent_dim = 'integer',
		sigma0 = 'numeric'
	),
	contains = 'vplot_model'
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
		membership = 'matrixOrNULL',
		centers = 'matrixOrNULL'
	),
	contains = 'vplot_autoencoder_model'
)

setClass(
	'vplot_autoencoder_cluster_v2_model',
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


