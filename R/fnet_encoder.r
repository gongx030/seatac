#' FnetEncoderModel
#' 
#' A FNet Encoder model for V-plot reconstruction
#'
#' @param latent_dim Latent dimension (default: 10L)
#' @param block_size Block size in base pairs (default: 640L)
#' @param bin_size Bin size in base pairs(default: 5L) 
#' @param filters0 Filter size after the latent layer (default: 128L)
#' @param fragment_size_range  Fragment size ranges (default: c(0L, 320L))
#' @param fragment_size_interval Fragment size interval (default: 10L)
#' @param rate Dropout rate (default: 0.1)
#' @param name Model name
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
FnetEncoderModel <- function(
	window_size = 2560L,
	patch_size = 80L,
	bin_size = 5L,
	embedding_dim = 256L,
	num_blocks = 4L,
	fragment_size_range  = c(0L, 320L),
	fragment_size_interval = 10L,
	mask_ratio = 0.9,
	rate = 0.1,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$bin_size <- bin_size
		self$n_bins_per_window <- as.integer(window_size / bin_size)
		self$n_bins_per_patch <- as.integer(patch_size / bin_size)

		self$fragment_size_range <- fragment_size_range
		self$fragment_size_interval <- fragment_size_interval
		br <- seq(fragment_size_range[1], fragment_size_range[2], by = fragment_size_interval)
		self$n_intervals <- length(br) - 1L
		self$breaks <- tf$constant(br)
		self$centers <- tf$constant((br[-1] + br[-length(br)]) / 2)
		self$positions <- tf$cast(seq(0 + bin_size / 2, window_size - bin_size / 2, by = bin_size) - (window_size / 2), tf$float32)
		self$embedding_dim <- embedding_dim
		self$num_patches <- as.integer(self$n_bins_per_window / self$n_bins_per_patch)
		self$n_patches_per_window <- as.integer(self$n_bins_per_window / self$n_bins_per_patch)

		self$mask <- MaskLayer(mask_ratio)

		self$patches_bin <- PatchLayer(
			patch_size = self$n_bins_per_patch,
			step_size = self$n_bins_per_patch
		)

		self$patches <- PatchLayer(
			patch_size = self$n_bins_per_patch * self$bin_size,
			step_size = self$n_bins_per_patch * self$bin_size
		)

		self$embed_vplots <- tf$keras$layers$Dense(embedding_dim)

#		self$embed_coverage <- tf$keras$layers$Dense(embedding_dim)

		self$positional_embedding <- tf$keras$layers$Embedding(input_dim = self$num_patches, output_dim = embedding_dim)
		
		self$encoder <- tf$keras$Sequential(lapply(1:num_blocks, function(i) FNetLayer(
			embedding_dim = embedding_dim,
			rate = rate
		)))

		self$output_vplot <- tf$keras$layers$Dense(self$n_intervals * self$n_bins_per_patch)
#		self$output_coverage <- tf$keras$layers$Dense(self$n_bins_per_patch * self$bin_size)

		function(x, ..., training = TRUE){

#			cvg <- x$coverage %>%
#				tf$expand_dims(1L) %>%
#				tf$expand_dims(3L) %>%
#				self$patches() %>%
#				self$embed_coverage()
					
			vplots <- x$vplots %>%
				self$patches_bin() %>%
				self$mask(training) %>%
				self$embed_vplots()

			pe <- tf$range(0L, self$num_patches) %>% 
				self$positional_embedding()

#			x <- vplots + pe + cvg
			x <- vplots + pe

			x <- x %>% self$encoder()

			vplots <- x %>% 
				self$output_vplot() %>%
				tf$reshape(shape(-1L, self$n_patches_per_window, self$n_intervals, self$n_bins_per_patch, 1L)) %>%
				tf$transpose(shape(0L, 2L, 1L, 3L, 4L)) %>%
				tf$reshape(shape(-1L, self$n_intervals, self$n_patches_per_window * self$n_bins_per_patch, 1L)) %>%
				tf$keras$activations$softmax(1L)

			list(
				vplots = vplots
			)
		}
	})
}
