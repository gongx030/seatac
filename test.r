library(roxygen2); library(devtools); devtools::create('analysis/seatac/packages/seatac')
library(roxygen2); library(devtools); devtools::document('analysis/seatac/packages/seatac')
    

# -----------------------------------------------------------------------------------
# [2019-05-31] Generate the fragment size distribution file for each ATAC-seq BAM file
# -----------------------------------------------------------------------------------
#tfe_enable_eager_execution(device_policy = 'silent')

library(tensorflow)
library(keras)
#use_backend(backend = "plaidml")

library(BSgenome.Mmusculus.UCSC.mm10)
filenames <- c(
  'ATAC_MEF_NoDox.bam', 
  'ATAC_MEF_Dox_D1.bam',
  'ATAC_MEF_Dox_D2.bam',
  'ATAC_MEF_Dox_D7.bam'
)
filenames <- sprintf('analysis/seatac/data/%s', filenames)
time_points <- factor(c('D0', 'D1', 'D2', 'D7'), c('D0', 'D1', 'D2', 'D7'))

# Etv2: chr7:30,604,535-30,664,933
#which <- GRanges(seqnames = 'chr7', range = IRanges(20000001, 40000000))
which <- GRanges(seqnames = 'chr7', range = IRanges(20000001, 80000000))
devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(filenames[1:2], which, genome = BSgenome.Mmusculus.UCSC.mm10, window_size = 10000, bin_size = 20, fragment_size_range = c(50, 500), fragment_size_interval = 25, epochs = 50, gpu = TRUE)

source('analysis/seatac/helper.r'); gr_file <- sprintf('%s/test.rds', PROJECT_DIR)
saveRDS(gr, file = gr_file)


par(mfrow = c(4, 1))
devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr[mcols(gr)$groups == 1], which = 'chr7:30,628,023-30,641,444')
plot(colMeans(mcols(gr)$predicted_counts), type = 'l')
plot(colMeans(mcols(gr)$counts), type = 'l')

devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr[mcols(gr)$groups == 1], which = 'chr7:30,638,023-30,641,444')


i <- 1

Xp <- model %>% predict(fs$X[i, , , drop = FALSE])
#image(Xp[1, , ], breaks = c(seq(0, 0.1, length.out = 100), 1), col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
image(Xp[1, , ], col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
y <- summary(as(fs$X[i, , ], 'dgCMatrix'))
points(y[, 1] / nrow(fs$X[i, , ]), y[, 2] / ncol(fs$X[i, , ]), pch = 3, cex = 1.25, col = 'black')


# -----------------------------------------------------------------------------------
# [2019-06-24] VAE
# -----------------------------------------------------------------------------------

library(tensorflow)
tfe_enable_eager_execution(device_policy = 'silent')
library(keras)
library(tfprobability)
library(tfdatasets)
library(futile.logger); flog.threshold(TRACE)
library(BSgenome.Mmusculus.UCSC.mm10)


# -----------------------------------------------------------------------------------
# [2019-07-11] segment the WT MEF
# -----------------------------------------------------------------------------------
#gs <- c('MEF_NoDox')
#gs <- c('MEF_Dox_D1')
#gs <- c('MEF_Dox_D2')
gs <- c('MEF_NoDox', 'MEF_Dox_D1')
window_size <- 320; bin_size <- 10
source('analysis/seatac/helper.r'); gr <- get_fragment_size(gs, window_size = window_size, bin_size = bin_size, min_reads_per_window = 10)

latent_dim <- 2; n_components <- 2; devtools::load_all('analysis/seatac/packages/seatac'); gr <- seatac(gr, latent_dim = latent_dim, n_components = n_components, prior = 'gmm', epochs = 50, batch_effect = FALSE, beta = 10)

source('analysis/seatac/helper.r'); gr_file <- sprintf('%s/results/test.rds', PROJECT_DIR)
saveRDS(gr, file = gr_file)


# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize the V-plot clusters for HMM prior
# -----------------------------------------------------------------------------------
library(gplots)
#par(mfrow = c(10, 1), mar = c(0.25, 2, 0.25, 2))
par(mfrow = c(5, 8), mar = c(2, 0.25, 2, 0.25))
X <- mcols(gr)$counts %>% as.matrix()
dim(X) <- c(length(gr), metadata(gr)$n_bins_per_window, metadata(gr)$n_intervals, metadata(gr)$num_samples)
X <- aperm(X, c(1, 4, 2, 3))
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window * metadata(gr)$n_intervals)
w <- 1 / rowSums(X); w[is.infinite(w)] <- 0
X <- as.matrix(Diagonal(x = w) %*% X)
dim(X) <- c(length(gr) * metadata(gr)$num_samples, metadata(gr)$n_bins_per_window,  metadata(gr)$n_intervals)
cls <- c(mcols(gr)$cluster)

lapply(1:metadata(gr)$model$n_components, function(k){
	X <- colSums(X[cls == k, , ], dims = 1)
	image(X, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE, main = k)
})

table(mcols(gr)$cluster[, 1], mcols(gr)$cluster[, 2])
table(mcols(gr)$cluster[, 2], mcols(gr)$cluster[, 3])
table(mcols(gr)$cluster[, 3], mcols(gr)$cluster[, 4])



# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize the V-plot clusters
# -----------------------------------------------------------------------------------
library(gplots)
#par(mfrow = c(10, 1), mar = c(0.25, 2, 0.25, 2))
par(mfrow = c(3, 5), mar = c(2, 0.25, 2, 0.25))
lapply(1:metadata(gr)$model$n_components, function(k) {
	X <- Reduce('+', lapply(1:metadata(gr)$num_samples, function(s){
		start <- 32 * 32 * (s - 1) + 1
		end <- 32 * 32 * s
		colSums(mcols(gr)$counts[mcols(gr)$cluster[, s] == k, start:end])
	}))
	X <- matrix(X, 32, 32)
	image(X, col = colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
})
table(mcols(gr)$cluster)


# -----------------------------------------------------------------------------------
# [2019-07-11] Compare the identified nucleosome and MNase-seq
# Need to run on the lab queue
# -----------------------------------------------------------------------------------
gs <- c('MEF_NoDox');  window_size <- 320; bin_size <- 10

peak_set <- sprintf('seatac_%s', paste(gs, collapse = '+'))
source('analysis/seatac/helper.r'); gr <- get_fragment_size(gs, window_size = window_size, bin_size = bin_size, min_reads_per_window = 10)
source('analysis/etv2_pioneer/helper.r'); peaks <- reCenterPeaks(gr, width = 1)
source('analysis/etv2_pioneer/helper.r'); bw_files <- get_bigwig_files()
gs <- c('MEF_NoDox_ATAC', 'MEF_nucleosome', 'MEF_MNase', 'MEF_H3', 'MEF_MNase2')
extend <- 320; w <- 10
source('analysis/etv2_pioneer/helper.r'); n2m_files <- normalizeToMatrix_batch(peaks, peak_set = peak_set, bw_files[gs], extend = extend, w = w, mc.cores = 4)

mat <- lapply(n2m_files, readRDS); names(mat) <- names(n2m_files)
cols <- c(
	'MEF_NoDox_ATAC' = 'red',
	'MEF_nucleosome' = 'blue',
	'MEF_MNase' = 'blue',
	'MEF_MNase2' = 'blue',
	'MEF_H3' = 'blue'
)
col_fun <- lapply(names(mat), function(g) colorRamp2(quantile(mat[[g]], c(0, 0.98)), c('white', cols[g])))
names(col_fun) <- names(n2m_files)

i <- sample.int(nrow(mat[[1]]), 5000)
#i <- 1:nrow(mat[[1]])
split_by <- mcols(peaks)$cluster[i]

axis_name <- c('-1k', 'summit', '+1k')
h <- EnrichedHeatmap(mat[['MEF_NoDox_ATAC']][i, ], split = split_by, col = col_fun[['MEF_NoDox_ATAC']], name = 'MEF_NoDox_ATAC', axis_name = axis_name, pos_line = FALSE) +
#EnrichedHeatmap(mat[['MEF_D7_ATAC']][i, ], col = col_fun[['MEF_D7_ATAC']], name = 'MEF_D7_ATAC', axis_name = axis_name, pos_line = FALSE)  +
#EnrichedHeatmap(mat[['EB_Dox_D25']][i, ], col = col_fun[['EB_Dox_D25']], name = 'EB_Dox_D25', axis_name = axis_name, pos_line = FALSE)  +
EnrichedHeatmap(mat[['MEF_nucleosome']][i, ], col = col_fun[['MEF_nucleosome']], name = 'MEF_nucleosome', axis_name = axis_name, pos_line = FALSE)  +
EnrichedHeatmap(mat[['MEF_MNase']][i, ], col = col_fun[['MEF_MNase']], name = 'MEF_MNase', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['MEF_MNase2']][i, ], col = col_fun[['MEF_MNase2']], name = 'MEF_MNase2', axis_name = axis_name, pos_line = FALSE)  + 
EnrichedHeatmap(mat[['MEF_H3']][i, ], col = col_fun[['MEF_H3']], name = 'MEF_H3', axis_name = axis_name, pos_line = FALSE) 

draw(h, heatmap_legend_side = 'right')



M <- as.matrix(mat[['MEF_NoDox_ATAC']])
plot(colMeans(M[mcols(peaks)$cluster == 1, ]), type = 'l', lwd = 2)
abline(v = 20.5, lwd = 2, lty = 2)


  mcols(win)$model$predicted_pattern <- lapply(1:n_components, function(k) colSums(Xp[cls == k, , , drop = FALSE], dims = 1))
  mcols(win)$model$observed_pattern <- V


n_components <- 10

df <- data.frame(cluster = factor(1:n_components), n = colSums(mcols(gr)$cluster))
ggplot(data = df, aes(x = cluster, y = n)) + geom_bar(stat = 'identity') + scale_y_log10()

lapply(1:10, function(k) image(metadata(gr)$observed_pattern[[k]], main = k))


# -----------------------------------------------------------------------------------
# [2019-06-26] Visualize a track that shows the following information:
# * Gene tracks
# * ATAC-seq coverage
# * fitted V-plot as well as observed counts
# * The cluster distribution
# -----------------------------------------------------------------------------------
library(gridExtra)
library(ggbio)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(SummarizedExperiment)
#devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr, which = 'chr7:30,628,023-30,641,444', txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)
#tx <- subsetByOverlaps(transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene), gr)

win <- Reduce('c', slidingWindows(GRanges(seqnames = 'chr7', range = IRanges(1, 145441459)), width = 1000, step = 1000))
A <- as.matrix(findOverlaps(win, gr))
A <- sparseMatrix(i = A[, 1], j = A[, 2], dims = c(length(win), length(gr)))
which <- win[order(rowSums(A), decreasing = TRUE)[50]]

devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr, which = resize(which, width = 10000, fix = 'center'), txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)




i <- 1:500 + 1000
image(as.matrix(mcols(gr)$counts[i, ]), col = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow'))
image(as.matrix(mcols(gr)$predicted_counts[i, ]), col = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow'))

# look at the clusters
lapply(1:15, function(k) image(colSums(mcols(gr)$predicted_counts[max.col(mcols(gr)$posterior) == k, , ], dims = 1), main = k))



devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr, which = 'chr7:30,628,023-30,641,444', txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)


which2 <- GRanges(seqnames = 'chr7', range = IRanges(10000001, 20000000))
devtools::load_all('analysis/seatac/packages/seatac'); gr <- model %>% predict(filenames[1], which, block_size = 10000000, min_reads_per_window = 20)



devtools::load_all('analysis/seatac/packages/seatac'); saveModel(model, dir = 'analysis/seatac/models/gmvae_20190625a')
devtools::load_all('analysis/seatac/packages/seatac'); model <- loadModel(dir = 'analysis/seatac/models/gmvae_20190625a')

which <- GRanges(seqnames = 'chr7', range = IRanges(10000001, 20000000))


		  lapply(1:n_components, function(k) image(colSums(fs$X[max.col(P) == k, , ], dims = 1), main = k))

			  i <- 1
			  par(mfrow = c(4, 6)); lapply(c(1, 50, 100, 200, 300, 500, 1000), function(i){image(fs$X[i, , ]); image(Xp[i, , ])})



par(mfrow = c(4, 1))
plot(colMeans(mcols(gr)$predicted_counts), type = 'l')
plot(colMeans(mcols(gr)$counts), type = 'l')

devtools::load_all('analysis/seatac/packages/seatac'); vplot(gr[mcols(gr)$groups == 1], which = 'chr7:30,638,023-30,641,444')


i <- 1

Xp <- model %>% predict(fs$X[i, , , drop = FALSE])
#image(Xp[1, , ], breaks = c(seq(0, 0.1, length.out = 100), 1), col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
image(Xp[1, , ], col = gplots::colorpanel(100, low = 'blue', mid = 'white', high = 'red'), axes = FALSE)
y <- summary(as(fs$X[i, , ], 'dgCMatrix'))
points(y[, 1] / nrow(fs$X[i, , ]), y[, 2] / ncol(fs$X[i, , ]), pch = 3, cex = 1.25, col = 'black')


; image(fs$X[i, , ])

  i <- 4
  plot(colMeans(Xp[i, , ]), main = 'predicted'); plot(colMeans(fs$X[i, , ]), main = 'observed')
	  plot(colMeans(Xp[1, , ]), main = 'predicted'); points(colMeans(Xp[2, , ]), col = 'red'); points(colMeans(Xp[3, , ]), col = 'blue')


; image(fs$X[i, , ])

  i <- 4
  plot(colMeans(Xp[i, , ]), main = 'predicted'); plot(colMeans(fs$X[i, , ]), main = 'observed')
	  plot(colMeans(Xp[1, , ]), main = 'predicted'); points(colMeans(Xp[2, , ]), col = 'red'); points(colMeans(Xp[3, , ]), col = 'blue')


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
#  observation_distribution = tfd_independent(tfd_bernoulli(logits = Theta, dtype = tf$float32)),
  observation_distribution = tfd_independent(tfd_normal(loc = Theta, scale = Theta)),
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

	    cvg_i <- sum(coverage(x)[bins])  # total reads coverage at each bin
	    names(cvg_i) <- NULL
			    cvg <- c(cvg, cvg_i)



H <- 4L
num_states <- 2   
num_steps <- 7
initial_state_logits <- rep(0, num_states)
transition_distribution <- tfd_categorical(
  probs = matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE) %>%
    tf$cast(tf$float32)
)
#observation_distribution <- tfd_multivariate_normal_diag(loc = matrix(0, 2, 4), scale_diag = rep(1, 4)) 
# We can combine these distributions into a single week long
# hidden Markov model with:
Theta <- tf$exp(tf$Variable(matrix(0, num_states, H), dtype = tf$float32))
d <- tfd_hidden_markov_model(
  initial_distribution = tfd_categorical(logits = initial_state_logits),
  transition_distribution = transition_distribution,
#  observation_distribution = tfd_independent(tfd_normal(loc = Theta, scale = Theta)),
  observation_distribution =  tfd_multivariate_normal_diag(loc = Theta, scale_diag = rep(1, H)),
  num_steps = num_steps
)
d$log_prob(d$sample(3L) %>% tf$reshape(shape(7, 3, 1, 4, 1)))




