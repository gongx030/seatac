#' fit
#'
setMethod(
	'fit',
	signature(
		model = 'vplot_glove_model',
		x = 'Vplots'
	),
	function(
		model,
		x,
		learning_rate = 1e-3, 
		batch_size = 256L,
		epochs = 100L
	){

		flog.info(sprintf('batch size(batch_size): %d', batch_size))

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		flog.info(sprintf('optimizer: Adam(learning_rate=%.3e)', learning_rate))

		Y <- t(x$counts) %*% x$counts
		u <- irlba(log(Y + 1), nu = model@latent_dim)$u

		z <- Diagonal(x = 1 / rowSums(x$counts)) %*% x$counts %*% u

		y_umap <- umap(as.matrix(z))$layout
		plot(y_umap, pch = 21, bg = classes)

	}
)

