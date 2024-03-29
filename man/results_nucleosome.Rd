% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.r
\name{results_nucleosome}
\alias{results_nucleosome}
\title{results_nucleosome}
\usage{
results_nucleosome(
  model,
  x,
  contrast,
  fragment_size_threshold = 150L,
  batch_size = 16L
)
}
\arguments{
\item{model}{a trained VaeModel object}

\item{x}{a Vplots object}

\item{contrast}{this argument species a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change (simplest case)}

\item{fragment_size_threshold}{Fragment size threshold for nucleosome reads}

\item{batch_size}{batch size (default: 16L)}
}
\value{
a GRanges object
}
\description{
Testing the difference of nucleosomes between two samples
}
