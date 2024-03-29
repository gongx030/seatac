% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vae.r
\name{VaeModel}
\alias{VaeModel}
\title{VaeModel}
\usage{
VaeModel(
  x,
  latent_dim = 10L,
  filters0 = 128L,
  filters = 32L,
  kernel_size = 3L,
  downsample_layers = 4L,
  upsample_layers = 4L,
  strides = c(2L, 2L),
  momentum = 0.8,
  rate = 0.1,
  name = NULL
)
}
\arguments{
\item{x}{a Vplots or a VplotsList object}

\item{latent_dim}{Latent dimension (default: 10L)}

\item{filters0}{Filter size after the latent layer (default: 128L)}

\item{filters}{Initial filter size of convolution layer (default: 32L)}

\item{kernel_size}{Kernel size in convolution and deconvolution  layers (default: 3L)}

\item{downsample_layers}{Downsample layers (default: 4L)}

\item{upsample_layers}{Upsample layers (default: 4L)}

\item{strides}{Convolution strides}

\item{momentum}{Momentum in BatchNormalization layer (default: 0.8)}

\item{rate}{Dropout rate (default: 0.1)}

\item{name}{Model name}
}
\value{
a VaeModel object
}
\description{
Build a VAE model for V-plot of multiple ATAC-seq datasets. This model takes the stacked V-plots of the same genomic regions as the input.
}
