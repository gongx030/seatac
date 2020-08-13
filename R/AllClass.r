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
	'vplot_parametric_vae_v3_vampprior_model',
	slot = c(
		num_pseudo_inputs = 'integer'
	),
	contains = 'vplot_parametric_vae_v3_model'
)


setClass(
	'vplot_parametric_vae_v3_gmm_model',
	slot = c(
		n_components = 'integer'
	),
	contains = 'vplot_parametric_vae_v3_model'
)

setClass(
	'vplot_parametric_vae_v4_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		mixture = 'kerastools.model.RModel',
		latent_dim = 'integer',
		sigma0 = 'numeric'
	),
	contains = 'vplot_model'
)


setClass(
	'vplot_parametric_vae_v5_model',
	slot = c(
		encoder = 'kerastools.model.RModel',
		prior = 'kerastools.model.RModel',
		decoder = 'kerastools.model.RModel',
		latent_dim = 'integer',
		mu = 'numeric',
		sigma = 'numeric',
		shape = 'numeric',
		scale = 'numeric'
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
	'vplot_autoencoder_v2_model',
	contains = 'vplot_autoencoder_model'
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

setClass(
	'vplot_glove_model',
	slot = c(
		latent_dim = 'integer'
	),
	contains = 'vplot_model'
)

setClass(
	'vplot_autoencoder_3d_model',
	slot = c(
		block_size = 'integer',
		n_bins_per_block = 'integer',
		n_blocks_per_window = 'integer'
	),
	contains = 'vplot_autoencoder_model'
)


setOldClass('rpytools.call.Encoder')
setOldClass('python.builtin.Encoder')
setClassUnion('rpytools.call.EncoderOrpython.builtin.Encoder', members = c('rpytools.call.Encoder', 'python.builtin.Encoder'))

setOldClass('rpytools.call.Decoder')
setOldClass('python.builtin.Decoder')
setClassUnion('rpytools.call.DecoderOrpython.builtin.Decoder', members = c('rpytools.call.Decoder', 'python.builtin.Decoder'))
setClass(
	'vplot_cvae_model',
	slot = c(
		encoder = 'rpytools.call.EncoderOrpython.builtin.Encoder',
		decoder = 'rpytools.call.DecoderOrpython.builtin.Decoder'
	)
)

