setOldClass('rpytools.call.Vplot2NucleosomeModel')
setOldClass('python.builtin.Vplot2NucleosomeModel')
setClassUnion('Vplot2NucleosomeModel', members = c('rpytools.call.Vplot2NucleosomeModel', 'python.builtin.Vplot2NucleosomeModel'))

Vplot2NucleosomeModel <- reticulate::PyClass(
	'Vplot2NucleosomeModel',
	inherit = tf$keras$Model,
	list(
		`__init__` = function(
			self,
			bin_size,
			latent_dim = 10L,
			block_size
		){
			super()$`__init__`()

			self$vplot_encoder <- VplotEncoder()
			self$dense_1 <- tf$keras$layers$Dense(units = latent_dim, activation = 'relu')
			self$dropout_1 <- tf$keras$layers$Dropout(rate = 0.1)

			self$dense_1 <- tf$keras$layers$Dense(units = bin_size, activation = 'relu')

			self$flatten <- tf$keras$layers$Flatten()

			NULL
		},
		call = function(self, x, training = TRUE){

			y <- x %>% 
				self$vplot_encoder() %>%
				self
			browser()
		}
	)
)

setMethod(
	'fit',
	signature(
		model = 'Vplot2NucleosomeModel',
		x = 'Vplots'
	),
	function(
		model,
		x,
		batch_size = 32L,
		epochs = 100L,
		min_reads_per_block = 10,
		test_size = 0.15
	){

		optimizer <- tf$keras$optimizers$Adam(1e-4, beta_1 = 0.9, beta_2 = 0.98, epsilon = 1e-9)
		mse <- tf$keras$losses$MeanSquaredError(reduction = 'none')

		# updating step
		train_step <- function(x, y){

			with(tf$GradientTape(persistent = TRUE) %as% tape, {
				y_pred <- model(x)
				loss <- mse(y, y_pred)  %>%
					tf$reduce_mean()
			})

			gradients <- tape$gradient(loss, model$trainable_variables)
			list(gradients, model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()

			loss

		} # train_step

		train_step <- tf_function(train_step) # convert to graph mode

		dataset <- select_blocks(x, 
			block_size = block_size, 
			min_reads = min_reads_per_block, 
			with_kmers = FALSE, 
			with_nucleoatac = TRUE, min_nucleoatac_signal = 0.01, max_nucleoatac_signal = 0.75
		) %>% 
			prepare_dataset(
				test_size = test_size, 
				batch_size = batch_size
			)

		for (epoch in seq_len(epochs)){

			loss_train <- 0
			loss_test <- 0

			iter <- make_iterator_one_shot(train_dataset)
			i <- 0
			until_out_of_range({
				batch <- iterator_get_next(iter)
				xi <- batch$vplot
				yi <- batch$nucleoatac
				loss_train <- loss_train + train_step(xi, yi)
				i <- i + 1
			})
			loss_train <- loss_train / i

			iter_test <- make_iterator_one_shot(test_dataset)
			i <- 0
			until_out_of_range({
				batch_test <- iterator_get_next(iter_test)
				xi <- batch_test$vplot
				yi <- batch_test$nucleoatac
				yi_pred <- model(xi)

				if (i == 0){
					tf$reduce_sum(yi, 0L) %>% plot(main = epoch)
					tf$reduce_sum(yi_pred, 0L)  %>% plot()
					xi %>% tf$reduce_sum(0L) %>% tf$squeeze(2L) %>% as.matrix() %>% t() %>% image()
				}

				loss_test <- loss_test + tf$reduce_mean(mse(yi, yi_pred))
				i <- i + 1
			})
			loss_test <- loss_test / i

			flog.info(sprintf('epoch=%6.d/%6.d | train_loss=%13.7f | test_loss=%13.7f', epoch, epochs, loss_train, loss_test))
		}

		model
	}
)
