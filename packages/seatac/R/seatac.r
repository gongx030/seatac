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
#' @import tfdatasets
NULL

#' seatac 
#'
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#'
#' @export
#'
#' @author Wuming Gong
#'
seatac <- function(filenames, which = NULL, genome, latent_dim = 10, num_states = 3, width = 200, expand = 2000, fragment_size_range = c(70, 500), fragment_size_interval = 10, min_reads_per_bin = 20, epochs = 5, batch_size = 256, sigma = 1){

  if (missing(filenames))
    stop('filenames are missing')

  num_samples <- length(filenames)

  existed <- file.exists(filenames)
  if (any(!existed))
    stop(sprintf('%s do not exist', paste(filenames[!existed], collapse = ', ')))

  is_pe <- sapply(filenames, testPairedEndBam)
  if(any(!is_pe)){
    stop(paste(filenames[!s_pe], collapse = ', '),"are not paired-end files.")
  }

  # genomic ranges covered by the BAM files
  gr <- Reduce('intersect', lapply(filenames, function(f){
    x <- idxstatsBam(f)
    GRanges(seqnames = x[, 'seqnames'], range = IRanges(1, x[, 'seqlength']))
  }))
  seqlevels(which, pruning.mode = 'coarse') <- seqlevels(gr)

  bins <- tileGenome(seqlengths(genome), tilewidth = width, cut.last.tile.in.chrom = TRUE)
  bins <- subsetByOverlaps(bins, which)
  seqlevels(bins, pruning.mode = 'coarse') <- seqlevels(gr)
  flog.info(sprintf('setting up %d overlapping bins (width=%d, expand=%d)', length(bins), width, expand))

  bins <- bins[1:10000]

  fs <- getFragmentSizeMatrix(filenames, bins, expand, fragment_size_range, min_reads_per_bin)
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

