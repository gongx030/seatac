library(roxygen2); library(devtools); devtools::create('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('analysis/seatac/packages/seatac')
    

# -----------------------------------------------------------------------------------
# [2019-05-31] Generate the fragment size distribution file for each ATAC-seq BAM file
# -----------------------------------------------------------------------------------
#tfe_enable_eager_execution(device_policy = 'silent')

library(BSgenome.Mmusculus.UCSC.mm10)
filenames <- c(
  'ATAC_MEF_NoDox.bam', 
  'ATAC_MEF_Dox_D1.bam',
  'ATAC_MEF_Dox_D2.bam',
  'ATAC_MEF_Dox_D7.bam'
)
filenames <- sprintf('analysis/seatac/data/%s', filenames)
time_points <- factor(c('D0', 'D1', 'D2', 'D7'), c('D0', 'D1', 'D2', 'D7'))

which <- GRanges(seqnames = 'chr7', range = IRanges(10000001, 20000000))
devtools::load_all('analysis/seatac/packages/seatac'); se <- seatac(filenames[1:2], which, genome = BSgenome.Mmusculus.UCSC.mm10, latent_dim = 20, window_size = 5000, bin_size = 10, epochs = 10)

sess <- tf$Session()
sess$run(tf$global_variables_initializer())
xx <- sess$run(tfb_softmax_centered()$forward(U$sample(1L)))

u <- sess$run(U$sample(1L))
x <- sess$run(tfb_softmax_centered(U)$sample(1L))


nc <- max(rowData(se)$cluster)
cols <- rainbow(nc)
cols2 <- col2rgb(cols)
cols2 <- sprintf('%d,%d,%d', cols2[1, ], cols2[2, ], cols2[3, ])

for (i in 1:2){
  se2 <- se[rowData(se)$group == i]
  f <- sprintf('~/Desktop/seatac_%d.bed', i)
  bed <- data.frame(
    chrom = seqnames(se2),
    chromStart = start(se2),
    chromEnd = end(se2),
    name = 1:length(se2),
    score = 1000,
    strand = '+',
    thickStart = start(se2),
    thickEnd = end(se2),
    itemRgb = cols2[rowData(se2)$cluster]
   )
  write.table(bed, f, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}


plot(colMeans(assays(se)$X[rowData(se)$cluster == 1, ]), type = 'l', col = cols[1], lwd = 1.5, ylim = c(0, 1))
lapply(2:4, function(k) lines(colMeans(assays(se)$X[rowData(se)$cluster == k, ]), type = 'l', col = cols[k], lwd = 1.5))

plot(colMeans(assays(se)$X_pred[rowData(se)$cluster == 1, ]), type = 'l', col = cols[1], lwd = 1.5, ylim = c(0, 1))
lapply(2:4, function(k) lines(colMeans(assays(se)$X_pred[rowData(se)$cluster == k, ]), type = 'l', col = cols[k], lwd = 1.5))


M <- ncol(assays(se)$X)

i <- 2
x <- hist(which(assays(se)$X[i, ] > 0), plot = FALSE)
xpd <- assays(se)$X_pred[i, ]; xpd <- xpd / sum(xpd)
plot(NA, xlim = c(0, 500), ylim = range(xpd))
lines(1:M, xpd, pch = 2, lwd = 2, type = 'l')
points(x$mids, x$density, lwd = 2, pch = 21, bg = 'gray', cex = 2)


  par(mfrow = c(4, 6))
  lapply(1:13, function(i)  plot(colMeans(Xp[igraph::membership(lou) == i, ]), main = i, type = 'l'))
#  lapply(1:num_states, function(i)  plot(colMeans(fs$X[igraph::membership(lou) == i, ]), main = i, type = 'l'))
  lapply(1:num_states, function(i)  plot(colMeans(Xp[cls == i, ]), main = i, type = 'l'))
  lapply(1:num_states, function(i)  plot(colMeans(fs$X[cls == i, ]), main = i, type = 'l'))

  flog.info(sprintf('clustering %d segments', nrow(Z)))
  k <- 100
  nn <- knn.index(Z, k)
  G <- as(sparseMatrix(i = rep(1:nrow(Z), k), j = c(nn), dims = c(nrow(Z), nrow(Z))), 'dgCMatrix')
  G <- graph.adjacency(G, mode = 'undirected')
  lou <- cluster_louvain(G)

  browser()


  Zp <- encoder %>% predict(X)

  Xp <- vae %>% predict(X)
  browser()
  browser()

devtools::load_all('analysis/seatac/packages/seatac'); se <- readFragmentSize(filenames[1], which = GRanges('chr1', range = IRanges(1, 1e7)), genome = BSgenome.Mmusculus.UCSC.mm10)

  tiles <- tileGenome(seqlengths(genome)['chr7'], tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
  tiles <- subsetByOverlaps(tiles, gr)
  seqlevels(tiles, pruning.mode = 'coarse') <- seqlevels(gr)
  flog.info(sprintf('dividing genome into %d non-overlapping tiles (tilewidth=%d)', length(tiles), tilewidth))
#  tiles <- tiles[1:10]

# -----------------------------------------------------------------------------------
# [2019-05-31] Testing Mixture model on simulated data
# -----------------------------------------------------------------------------------
library(tensorflow)
library(tfprobability)
library(gtools)
library(futile.logger)
library(Matrix)
#tfe_enable_eager_execution()

H <- 10L # number of segments
K <- 3L # nubmer of clusters
M <- 200  # number of samples
set.seed(1)
P0 <- rdirichlet(K, rep(0.1, H))  
N <- sample(20:100, M, replace = TRUE)  # number of total reads in each segment
cluster <- sample(1:K, M, replace = TRUE)
X <- t(sapply(1:M, function(m) rmultinom(1, N[m], P0[cluster[m], ])))
X <- as.matrix(Diagonal(x = 1 / rowSums(X)) %*% X)


K <- 3L  # number of clusters
cat_logits <- tf$Variable(rep(0, K), dtype = tf$float32)
Theta <- tf$Variable(matrix(log(colMeans(X) + 1), K, H, byrow = TRUE) + rnorm(K * H), dtype = tf$float32)
components_distribution <- tfd_independent(tfd_bernoulli(logits = Theta, dtype = tf$float32))
mix <- tfd_mixture_same_family(
  mixture_distribution = tfd_categorical(logits = cat_logits),
  components_distribution = components_distribution
)
total_log_prob <- mix %>% tfd_log_prob(as.matrix(X))
total_log_prob <- tf$reduce_sum(tf$multiply(total_log_prob, N))

train_op <- tf$train$AdamOptimizer(0.5)$minimize(-total_log_prob)

for (i in 1:1000){
  sess$run(train_op)
  if (i %% 100 == 0){
    flog.info(sprintf('step=%d | log prob=%.3f', i, sess$run(total_log_prob)))
    Theta_p <- sess$run(Theta)
    print(Theta_p)
  }
}

x_post <- tfd_independent(tfd_bernoulli(components_distribution$mean()))
post <- sess$run(x_post$log_prob(tf$reshape(X, shape(M, 1, H))))


# -----------------------------------------------------------------------------------
# [2019-05-31] Testing Mixture model on simulated data
# -----------------------------------------------------------------------------------
library(tensorflow)
library(tfprobability)
library(gtools)
library(futile.logger)
library(Matrix)
#tfe_enable_eager_execution()

H <- 10L # number of segments
K <- 3L # nubmer of clusters
M <- 200  # number of samples


# HMM
H <- 10L
p0 <- 0.05
num_states <- 3
K <- 2
num_steps <- 100
transition_probs <- matrix(p0 / (num_states - 1), num_states, num_states)
diag(transition_probs) <- 1 - p0
Lambda <- rbind(c(1, 2, 3, 4, 5), c(5, 4, 3, 2, 1), c(1, 1, 1, 1, 1))
X <- matrix(sample(c(0, 1), num_steps * H, replace = TRUE), num_steps, H)
Theta <- tf$exp(tf$Variable(matrix(log(colMeans(X) + 1), num_states, H, byrow = TRUE) + rnorm(num_states * H), dtype = tf$float32))
transition_distribution <- tfd_categorical(probs = transition_probs %>% tf$cast(tf$float32))
initial_state_logits <- rep(0, num_states)
hmm <- tfd_hidden_markov_model(
  initial_distribution = tfd_categorical(logits = initial_state_logits),
  transition_distribution = transition_distribution,
  observation_distribution = tfd_independent(tfd_bernoulli(logits = Theta, dtype = tf$float32)),
  num_steps = num_steps
)

Y2 <- tf$reshape(hmm$sample(2L), shape(100, 2, 1, 10, 1))
hmm$log_prob(Y2)
hmm$posterior_marginals(Y2)

hmm$log_prob(X)
hmm %>% tfd_log_prob(X)





      browser()


total_log_prob <- tf$reduce_sum(%>% tfd_log_prob(Y))

num_states <- 3
num_steps <- 200L

total_log_prob <- 
state <- 1
for (i in 2:num_steps) state[i] <- sample.int(num_states, 1, prob = transition_probs[state[i - 1], ])


p0 <- 0.05
transition_probs <- matrix(p0 / (num_states - 1), num_states, num_states)
diag(transition_probs) <- 1 - p0
Lambda <- rbind(c(1, 2, 3, 4, 5), c(5, 4, 3, 2, 1), c(1, 1, 1, 1, 1))
Theta <- tf$exp(tf$Variable(matrix(log(colMeans(X) + 1), K, H, byrow = TRUE) + rnorm(K * H), dtype = tf$float32))
transition_distribution <- tfd_categorical(probs = transition_probs %>% tf$cast(tf$float32))
initial_state_logits <- rep(0, num_states)
hmm <- tfd_hidden_markov_model(
  initial_distribution = tfd_categorical(logits = initial_state_logits),
  transition_distribution = transition_distribution,
  observation_distribution = tfd_independent(tfd_bernoulli(logits = Theta, dtype = tf$float32),
#  observation_distribution = tfd_independent(tfd_poisson(rate = Theta)),
#  observation_distribution = tfd_poisson(rate = Theta[, 1]),
#  observation_distribution = tfd_multinomial(10L, logits = rbind(1:4, 5:8, 9:12) %>% tf$cast(tf$float32)),
  num_steps = num_steps
)



Y <- X %>% tf$cast(tf$float32) %>% tf$reshape(shape(num_steps, M))
total_log_prob <- hmm %>% tfd_log_prob(Y)

Y <- X %>% tf$cast(tf$float32) %>% tf$reshape(shape(M, num_steps, 1))
posterior_dists <- hmm$posterior_marginals(Y)
sess <- tf$Session()
posterior_probs_ <- sess$run(posterior_dists$probs)  # extract the probabilities.



hmm$posterior_mode(Y)


  flog.info('initialize the HMM model')
  p0 <- 0.05
  transition_probs <- matrix(p0 / (num_states - 1), num_states, num_states)
  diag(transition_probs) <- 1 - p0
  transition_distribution <- tfd_categorical(probs = transition_probs %>% tf$cast(tf$float32))
  initial_state_logits <- rep(0, num_states)

#  Mu <- tf$Variable(matrix(0, num_states, M, byrow = TRUE) + rnorm(M * num_states), dtype = tf$float32)
#  Sigma <- tf$Variable(matrix(1, num_states, M, byrow = TRUE), dtype = tf$float32)

  rate_prior <- tfd_log_normal(1, 1)

  # simulate data
  Lambda_sim <- log(rbind(1:M, M:1))
  components_distribution <- tfd_independent(tfd_poisson(log_rate = Lambda_sim))
  hmm <- tfd_hidden_markov_model(
    initial_distribution = tfd_categorical(logits = initial_state_logits),
    transition_distribution = transition_distribution,
    observation_distribution = components_distribution,
    num_steps = n_stages
  )
  sess <- tf$Session()
  sess$run(tf$global_variables_initializer())
  X_sim <- sess$run(hmm$sample(250L))
  X <- X_sim %>% tf$cast(tf$float32)

  #posterior_dists <- hmm$posterior_marginals(X)
  #posterior_probs_ <- sess$run(hmm$posterior_marginals(X)$probs)  # extract the probabilities.

  Lambda <- tf$exp(tf$Variable(matrix(0, num_states, M) + rnorm(num_states * M), dtype = tf$float32))
  components_distribution <- tfd_independent(tfd_poisson(rate = Lambda))
  hmm <- tfd_hidden_markov_model(
    initial_distribution = tfd_categorical(logits = initial_state_logits),
    transition_distribution = transition_distribution,
    observation_distribution = components_distribution,
    num_steps = n_stages
  )
  latent_prior <- tfd_multivariate_normal_diag(loc  = tf$zeros(latent_dim), scale_identity_multiplier = 1)
  latent_prior <- tfd_independent(tfd_normal(loc  = tf$zeros(latent_dim), scale = 1))
    
  input_layer <- layer_input(shape = M)
  encoder <- input_layer %>% 
    layer_dense(units = params_size_multivariate_normal_tri_l(latent_dim)) %>%
    layer_multivariate_normal_tri_l(event_size = latent_dim) %>%
    layer_kl_divergence_add_loss(
      distribution = tfd_independent(
        tfd_normal(loc = rep(0, latent_dim), scale = rep(1, latent_dim)),
        reinterpreted_batch_ndims = 1
      ),
      weight = 1
    )
    

    keras_model_custom(name = name, function(self){
      self$dense <- layer_dense(units = latent_dim * 2)
      function(x, mask = NULL){
        x <- x %>% self$dense()
        tfd_multivariate_normal_diag(
          loc = x[, 1:latent_dim], 
          scale_diag = tf$nn$softplus((x[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)),
          activity_regularizer = tfp_
        )
      }
    })
  }

  encoder <- encoder_model()
  iter <- make_iterator_one_shot(train_dataset)
  x <-  iterator_get_next(iter)
  approx_posterior <- encoder(x)
  approx_posterior_sample <- approx_posterior$sample()

  decoder_model <- function(name = NULL){
    keras_model_custom(name = name, function(self){
      self$dense <- layer_dense(units = M, activation = 'relu')
      function(x, mask = NULL){
        x <- x %>% self$dense()
        x
      }
    })
  }
  decoder <- decoder_model()
  decoder_likelihood <- decoder(approx_posterior_sample)

  optimizer <- tf$train$AdamOptimizer(1e-4)
  compute_kl_loss <- function(
    latent_prior,
    approx_posterior,
    approx_posterior_sample
  ){
    kl_div <- approx_posterior$log_prob(approx_posterior_sample) - latent_prior$log_prob(approx_posterior_sample)
    avg_kl_div <- tf$reduce_mean(kl_div)
    avg_kl_div
  }




  browser()


  P <- tf$nn$softmax(tf$matmul(U, V, transpose_a = TRUE))
  Pi <- tfd_categorical(logits = matrix(0, N, num_states), name = 'Pi')
#  multi <- tfd_multinomial(total_count %>% tf$reshape(shape(N, 1)), probs = P)
  multi <- tfd_multinomial(total_count = , probs = P)
  X <- tfd_mixture_same_family(mixture_distribution = Pi, components_distribution = multi, name = 'X')


  browser()

  model$vae %>% fit(fs$X, fs$X, shuffle = TRUE, epochs = epochs, batch_size = batch_size, validation_split = 0.1)
  Z <- model$encoder %>% predict(fs$X, batch_size = batch_size)
  Xp <- model$vae %>% predict(fs$X, batch_size = batch_size)

  flog.info(sprintf('clustering %d segments', nrow(Z)))
  cls <- kmeans(Z, num_states)$cluster
  se <- SummarizedExperiment(assays = list(X = fs$X, X_pred = Xp), rowRanges = fs$bins)
  rowData(se)$cluster <- cls
  rowData(se)$group <- fs$group
  se

} # seatac


  U <- tf$Variable(matrix(rnorm(num_states * M), num_states, M), dtype = tf$float32)
  P <- tf$nn$softmax(U)
  Y <- tf$lgamma(total_count) - tf$reduce_sum(tf$lgamma(X + 1) - tf$multiply(X, tf$log(P + 1e-5)), axis = 2L)
  V <- tf$Variable(matrix(rnorm(N * num_states), N, num_states), dtype = tf$float32)
  Q <- tf$nn$softmax(V)
  loss <- - (tf$reduce_sum(Q * (Y + sigma * tf$log(Q + 1e-5))) + 100 * tf$norm(U))

  train <- tf$train$AdamOptimizer(0.5)$minimize(loss)

  sess <- tf$Session()
  sess$run(tf$global_variables_initializer())

  for (epoch in 1:epochs){
    sess$run(train)
    if (epoch %% 10 == 0){
      flog.info(sprintf('epoch=%d/%d | loss=%.3f', epoch, epochs, sess$run(loss)))
    }
  }

  P_pred <- sess$run(P)
  Q_pred <- sess$run(Q)
  lapply(1:num_states, function(k) plot(P_pred[k, ], type = 'l', main = k))
  lapply(1:num_states, function(k) plot(colMeans(fs$X[max.col(Q_pred) == k, ]), type = 'l', main = k))
  table(max.col(Q_pred))


  browser()
  Up <- sess$run(U)


  next_batch %>% tf$reshape(shape(NULL, 3, 22))



  multinom <- tfd_multinomial(total_count %>% tf$reshape(shape(N, 1)) %>% tf$cast(tf$float64), probs = P)
  multinom <- tfd_multinomial(10, probs = P)
  X <- tfd_mixture_same_family()

  M <- ncol(fs$X)
  N <- nrow(fs$X)
  train_dataset <- fs$X %>% 
    as.matrix() %>% 
    tf$cast(tf$float32) %>% 
    tensor_slices_dataset() %>% 
    dataset_batch(batch_size)

  encoder_model <- function(latent_dim, name = NULL){
    keras_model_custom(name = name, function(self){
      self$dense <- layer_dense(units = 2 * latent_dim)
      function(x, mask = NULL){
        x <- x %>% self$dense()
        tfd_multivariate_normal_diag(
          loc = x[, 1:latent_dim],
          scale_diag = tf$nn$softplus(x[, (latent_dim + 1):(2 * latent_dim)] + 1e-5)
        ) 
      }
    })
  }
  encoder <- encoder_model(latent_dim)

  decoder_model <- function(output_dim, name = NULL){
    keras_model_custom(name = name, function(self){
      self$dense <- layer_dense(units = output_dim)
      function(x, mask = NULL){
        x[[1]] <- self$dense(x[[1]])
        tfd_multinomial(total_count = x[[2]], logits = x[[1]])
      }
    })
  }
  decoder <- decoder_model(M)

  optimizer <- tf$train$AdamOptimizer(1e-4)
  compute_kl_loss <- function(latent_prior, approx_posterior, approx_posterior_sample) {
    kl_div <- approx_posterior$log_prob(approx_posterior_sample) - latent_prior$log_prob(approx_posterior_sample)
    avg_kl_div <- tf$reduce_mean(kl_div)
    avg_kl_div
  }

  latent_prior <- tfd_multivariate_normal_diag(loc = rep(0, latent_dim), scale_identity_multiplier = 1)

  for (epoch in seq_len(epochs)) {

    iter <- make_iterator_one_shot(train_dataset)

    total_loss <- 0
    total_loss_nll <- 0
    total_loss_kl <- 0

    until_out_of_range({

      x <-  iterator_get_next(iter)
      total_count <- tf$reduce_sum(x, axis = 1L)

      with(tf$GradientTape(persistent = TRUE) %as% tape, {

        approx_posterior <- encoder(x)
        approx_posterior_sample <- approx_posterior$sample()
        decoder_likelihood <- decoder(list(approx_posterior_sample, total_count))

        nll <- -decoder_likelihood$log_prob(x)
        avg_nll <- tf$reduce_mean(nll)

        kl_loss <- compute_kl_loss(latent_prior, approx_posterior, approx_posterior_sample)

        loss <- kl_loss + avg_nll
      })

      total_loss <- total_loss + loss
      total_loss_nll <- total_loss_nll + avg_nll
      total_loss_kl <- total_loss_kl + kl_loss

      encoder_gradients <- tape$gradient(loss, encoder$variables)
      decoder_gradients <- tape$gradient(loss, decoder$variables)

      optimizer$apply_gradients(
        purrr::transpose(list(encoder_gradients, encoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )

      optimizer$apply_gradients(
        purrr::transpose(list(decoder_gradients, decoder$variables)),
        global_step = tf$train$get_or_create_global_step()
      )
    })
    flog.info(sprintf('epoch=%d/%d | negative likelihood=%.3f | kl=%.3f | total=%.3f', epoch, epochs, total_loss_nll, total_loss_kl, total_loss))
  }
  browser()

  Z <- encoder(fs$X %>% as.matrix())$loc %>% as.matrix()
  cls <- kmeans(Z, num_states)$cluster
  lapply(1:num_states, function(k) plot(colMeans(fs$X[cls == k, ]), type = 'l', main = k))

}

  sampling <- function(arg){
    z_mean <- arg[, 1:(latent_dim)]
    z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]),
      mean = 0.,
      stddev = 1.0
    )
    z_mean + k_exp(z_log_var / 2) * epsilon
  }

cae <- function(input_dim, feature_dim, latent_dim = 10){

  input_layer <- layer_input(shape = shape(input_dim, feature_dim))
  encoder <- input_layer %>% bidirectional(layer_gru(units = latent_dim))

  decoder <- encoder %>%
    layer_repeat_vector(input_dim) %>%
    layer_gru(latent_dim, return_sequences = TRUE) %>%
    time_distributed(layer_dense(units = feature_dim, activation = 'softmax'))

  model <- keras_model(input_layer, decoder)
  model %>% compile(optimizer = 'rmsprop', loss = 'binary_crossentropy')
  model

} #  vae
