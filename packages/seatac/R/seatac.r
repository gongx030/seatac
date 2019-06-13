#' seatac
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom matrixStats rowSds rowVars rowMedians
#' @import tensorflow
#' @import tfprobability 
#' @import keras
#' @importFrom futile.logger flog.info
#' @importFrom GenomicRanges tileGenome resize intersect reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools testPairedEndBam ScanBamParam scanBamFlag idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom abind abind
NULL

#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(filenames, time_points, genome, latent_dims = 10, num_states = 2, width = 200, expand = 2000, tilewidth = 1e7, fragment_size_range = c(70, 750), fragment_size_interval = 10, min_reads_per_bin = 20, epochs = 5){

  if (missing(filenames))
    stop('filenames are missing')

  if (missing(time_points))
    stop('time_points are missing')

  num_samples <- length(filenames)
  num_stages <- length(time_points)

  # number of segments
  M <- (fragment_size_range[2] - fragment_size_range[1]) / fragment_size_interval + 1

  existed <- file.exists(filenames)
  if (any(!existed))
    stop(sprintf('%s do not exist', paste(filenames[!existed], collapse = ', ')))

  is_pe <- sapply(filenames, testPairedEndBam)
  if(any(!is_pe)){
    stop(paste(filenames[!s_pe], collapse = ', '),"are not paired-end files.")
  }

  if (!is.factor(time_points))
    stop('time_points must be factors')

  if (length(time_points) != length(filenames))
    stop('filenames and time_points must have the same length')

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
   }))

  tiles <- tileGenome(seqlengths(genome), tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
  tiles <- subsetByOverlaps(tiles, gr)
  seqlevels(tiles, pruning.mode = 'coarse') <- seqlevels(gr)
  flog.info(sprintf('dividing genome into %d non-overlapping tiles (tilewidth=%d)', length(tiles), tilewidth))
  tiles <- tiles[1]

  bins <- tileGenome(seqlengths(genome), tilewidth = width, cut.last.tile.in.chrom = TRUE)
  bins <- subsetByOverlaps(bins, tiles)
  seqlevels(bins, pruning.mode = 'coarse') <- seqlevels(gr)
  flog.info(sprintf('dividing genome into %d overlapping bins (width=%d, expand=%d)', length(bins), width, expand))

  sampling <- function(arg){
    z_mean <- arg[, 1:(latent_dims)]
    z_log_var <- arg[, (latent_dims+ 1):(2 * latent_dims)]
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]), 
      mean = 0.,
      stddev = 1.0
    )
    z_mean + k_exp(z_log_var / 2) * epsilon
  }

  get_gamma <- function(z){
  }

  pii <- tf$Variable(rep(1 / num_states, num_states), dtype = tf$float32)  # prior for membership
  U <- tf$Variable(matrix(0, latent_dims, num_states), dtype = tf$float32)  # cluster mean
  Lambda <- tf$Variable(matrix(1, latent_dims, num_states), dtype = tf$float32) # cluster variance

  for (epoch in 1:epochs){
    for (g in 1:length(tiles)){ # for each tile

      flog.info(sprintf('epoch=%d/%d | tile=%d/%d | %s:%d-%d', epoch, epochs, g, length(tiles), seqnames(tiles[g]), start(tiles[g]), end(tiles[g])))

      bins2 <- bins[bins %within% tiles[g]]
      bins2 <- resize(bins2, expand, fix = 'center')
      bins2 <- trim(bins2)

      X <- lapply(1:num_samples, function(i) readFragmentSize(filenames[i], bins2, fragment_size_range, fragment_size_interval))
      included <- rowSums(do.call('cbind', lapply(X, rowSums)) > min_reads_per_bin) == num_samples
      X <- lapply(1:num_samples, function(i) X[[i]][included, ])
      X <- abind(lapply(X, as.matrix), along = 1.5)  # combine a list of n_intervals ~ n_segment matrices into a n_intervals ~ n_stages ~ n_segmentsarray
      X <- X[, 1, ]

      input_layer <- layer_input(shape = shape(M))
      encoder_layers <- input_layer %>% 
        layer_reshape(target_shape = shape(M, 1)) %>% 
        layer_conv_1d(filters = 5, kernel_size = 10,  activation = 'relu', padding = 'same') %>%
        layer_max_pooling_1d(pool_size = 4) %>% 
        layer_flatten() 

      z_mean <- layer_dense(encoder_layers, latent_dims)
      z_log_var <- layer_dense(encoder_layers, latent_dims)
      z <- layer_concatenate(list(z_mean, z_log_var)) %>% layer_lambda(sampling)
      decoder_layers <- layer_dense(units = M, activation = 'relu')
      output_layer <- decoder_layers(z)

      browser()

      vae <- keras_model(input_layer, output_layer)

      vae_loss <- function(x, x_decoded_mean){

        xent_loss <- (M / 1.0) * loss_mean_squared_error(x, x_decoded_mean) # the reconstruction loss
        kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
        xent_loss + kl_loss
      }

      vae %>% compile(optimizer = 'rmsprop', loss = vae_loss)

      vae %>% fit(X, X, shuffle = TRUE, epochs = 10, batch_size = 256)
      Xp <- vae %>% predict(X)

      encoder <- keras_model(input_layer, z_mean)
      Zp <- encoder %>% predict(X)

      Xp <- vae %>% predict(X)
      
      browser()


      basis <- k_random_normal(shape = shape(num_states, M))
      basis <- k_exp(basis)
      Y <- k_dot(V, basis)

      model <- keras_model(list(input_layer), Y)
      model_optimizer <- optimizer_adam(lr = 1e-2)
      model %>% compile(optimizer = model_optimizer, loss = 'mse')



#      X <- log(X + 1)

#      N <- 10000L
#      X <- X[sample.int(dim(X)[1], N), , ]
      N <- dim(X)[1]

      #S <- tf$Variable(array(rnorm(N * num_stages * num_states), dim = c(N, num_stages, num_states)), dtype = tf$float32) # stage ~ cluster
#      P <- tf$nn$softmax(S)

      U <- tf$Variable(array(rnorm(N * num_stages * num_states), dim = c(N, num_stages, num_states)) %>% tf$cast(tf$float32))
      U <- tf$exp(U)
      V <- tf$Variable(array(rnorm(num_states * M), dim = c(num_states, M)) %>% tf$cast(tf$float32))
      V <- tf$exp(V)
      pois <- tfd_poisson(tf$matmul(U, V))
      J1 <- -pois$log_prob(X) %>% tf$reduce_sum()

      # the segment weight
#      b_interval <- tf$Variable(rnorm(N), dtype = tf$float32)
 #     b_interval <- b_interval %>% tf$reshape(shape(N, 1L, 1L, 1L))

  #    Xe <- b_interval + b_stage + Beta + Theta
       
  #    Xd <- X %>% tf$cast(tf$float32) %>% tf$reshape(shape(N, num_stages, 1L, M)) - Xe
  #    J1 <- Xd %>% tf$square() %>% tf$reduce_sum(axis = 3L) %>% tf$multiply(P) %>% tf$reduce_sum()
#      J2 <- -tf$reduce_sum(P[, 1, ] * P[, 2, ] + P[, 2, ] * P[, 3, ] + P[, 3, ] * P[, 4, ])
      loss <- J1

      train_op <- tf$train$AdamOptimizer(0.1)$minimize(loss)
      sess <- tf$Session()
      sess$run(tf$global_variables_initializer())

      for (i in 1:1500){
        sess$run(train_op)
        if (i %% 100 == 0){
#          flog.info(sprintf('step=%d | J1=%.3f | J2=%.3f', i, sess$run(J1), sess$run(J2)))
          flog.info(sprintf('step=%d | J1=%.3f | J2=%.3f', i, sess$run(J1), 0))
        }
      }
      browser()

      plot(sess$run(b_stage), sess$run(tf$reduce_sum(X, c(0L, 2L))))
      plot(sess$run(b_interval), sess$run(tf$reduce_sum(X, c(1L, 2L))))
      lapply(1:num_stages, function(stage){
        plot(sess$run(Beta)[stage, 1, ], type = 'l', main = sprintf('stage %d', stage))
      })
      lapply(1:num_states, function(state){
        plot(sess$run(Theta)[1, state, ], type = 'l', main = sprintf('state %d', state))
      })
      lapply(1:num_states, function(state){
        plot(sess$run(V[state, ]), type = 'l', main = sprintf('state %d', state))
      })

      clusters <- sess$run(tf$argmax(P, axis = 2L))
      lapply(1:num_stages, function(stage){
        plot(colMeans(X[clusters[, stage] == 0, stage, ]), type = 'l', main = sprintf('stage %d', stage), ylim = c(0, 2))
        lapply(2:num_states, function(state) lines(colMeans(X[clusters[, stage] == state - 1, stage, ]), col = state))
      })


      Beta_p <- sess$run(Beta)

      Theta_p <- sess$run(Theta)


      lapply(1:num_stages, function(stage){
        plot(colMeans(X[clusters[, 1] == 0, stage, ]), type = 'l', main = 'stage', ylim = range(X[, stage, ]))
        lines(colMeans(X[clusters[, 1] == 1, stage, ]), col = 'red')
      })

      browser()
      table(factor(clusters[, 1], 0:(num_states - 1)), factor(clusters[, 2], 0:(num_states - 1)))
      table(factor(clusters[, 2], 0:(num_states - 1)), factor(clusters[, 3], 0:(num_states - 1)))
      table(factor(clusters[, 3], 0:(num_states - 1)), factor(clusters[, 4], 0:(num_states - 1)))

#      N <- sapply(X, rowSums) # total number of reads in each segment
#      BG <- lapply(X, function(x) colSums(x) / sum(x))# the background distribution of fragment size
#      X <- lapply(1:n_samples, function(i) Diagonal(x = 1 / N[, i]) %*% X[[i]])  # scale to each segment vector to prob vector
#      # normalize the segment fragment size distribution against the background distribution
#      X <- lapply(1:n_samples, function(i) as.matrix(log( (X[[i]] + 1) %*% Diagonal(x = 1 / (BG[[i]] + 1))) ))

#      X0 <- X[1:250, , ]
#      X2 <- tf$constant(X0, dtype = tf$float32)

#      X2 <- hmm$sample(250L)
#      hmm$log_prob(tf$reshape(X2, shape(nrow(X2), n_stages, 1, M, 1)))
#      total_log_prob <- tf$reduce_sum(hmm %>% tfd_log_prob(X2))
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X2, shape(nrow(X2), M, 1, n_stages, 1))))

#     total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(nrow(X), n_stages, 1, M, 1)))) + tf$reduce_sum(rate_prior$log_prob(Lambda))
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[2], 1, dim(X)[3], 1)))) + tf$reduce_sum(rate_prior$log_prob(Lambda))
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[2], 1, dim(X)[3], 1)))) 
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[2], 1, dim(X)[3], 1)))) 
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[3], 1, dim(X)[2], 1)))) 
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[3], dim(X)[2], 1))))
#      total_log_prob <- hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[2], dim(X)[3], 1))) +  tf$reduce_sum(rate_prior$log_prob(Lambda))
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(dim(X)[1], dim(X)[2] * dim(X)[3], 1))))
#      total_log_prob <- tf$reduce_sum(hmm$log_prob(tf$reshape(X, shape(n_stages, dim(X)[1], M, 1))))


      Lambda_pred <- sess$run(Lambda)
      browser()


      browser()

      posterior_dists <- hmm$posterior_marginals(X)
      posterior_probs_ <- sess$run(hmm$posterior_marginals(X)$probs)  # extract the probabilities.
      posterior_probs_


      total_log_prob <- hmm$log_prob(rpois(100, 5))

#      hmm %>% tfd_log_prob

      y <- tf$placeholder(tf$float32, as.matrix(X[[1]][1:100, ]))
      total_log_prob <- hmm$log_prob(as.matrix(X[[1]][1:100, ]))




      

      tfd_mixture(cat = tfd_categorical(logits = initial_state_logits), components = lapply(Lambda, function(lambda) tfd_poisson(lambda)))




    }
  }
  browser()
}



