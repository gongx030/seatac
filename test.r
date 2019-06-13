library(roxygen2); library(devtools); devtools::create('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('analysis/seatac/packages/seatac')
    

# -----------------------------------------------------------------------------------
# [2019-05-31] Generate the fragment size distribution file for each ATAC-seq BAM file
# -----------------------------------------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
filenames <- c(
  'ATAC_MEF_NoDox.bam', 
  'ATAC_MEF_Dox_D1.bam',
  'ATAC_MEF_Dox_D2.bam',
  'ATAC_MEF_Dox_D7.bam'
)
filenames <- sprintf('analysis/seatac/data/%s', filenames)
se_files <- gsub('.bam', '.fragment_size.rds', filenames)
time_points <- factor(c('D0', 'D1', 'D2', 'D7'), c('D0', 'D1', 'D2', 'D7'))

devtools::load_all('analysis/seatac/packages/seatac'); seatac(filenames, time_points, genome = BSgenome.Mmusculus.UCSC.mm10, latent_dims = 10, num_states = 3)

devtools::load_all('analysis/seatac/packages/seatac'); se <- readFragmentSize(filenames[1], which = GRanges('chr1', range = IRanges(1, 1e7)), genome = BSgenome.Mmusculus.UCSC.mm10)


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
sess <- tf$Session()
sess$run(tf$global_variables_initializer())

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
