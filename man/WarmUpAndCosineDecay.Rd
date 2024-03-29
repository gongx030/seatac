% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scheduler.r
\name{WarmUpAndCosineDecay}
\alias{WarmUpAndCosineDecay}
\title{Warmup and cosine decay scheduler, the schedular used in the SimCLR paper}
\usage{
WarmUpAndCosineDecay(
  base_learning_rate = 0.3,
  num_examples = NULL,
  batch_size = 256L,
  warmup_epochs = 10,
  epochs = 100L
)
}
\arguments{
\item{base_learning_rate}{Baseline learning rate (default: 0.3)}

\item{num_examples}{Number of total training samples}

\item{batch_size}{Batch size (default : 256)}

\item{warmup_epochs}{Warmup epochs (default: 10)}

\item{epochs}{Total training epochs}
}
\value{
a WarmUpAndCosineDecay object
}
\description{
Adopted from https://github.com/google-research/simclr/blob/67b562524223e65c0414a15c3bd7ac007d28f629/tf2/model.py#L78
}
