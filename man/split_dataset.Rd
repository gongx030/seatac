% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{split_dataset}
\alias{split_dataset}
\title{split_dataset}
\usage{
split_dataset(x, test_size = 0.15, batch_size = 64L)
}
\arguments{
\item{x}{a tfdataset object}

\item{test_size}{The ratio of the testing set (default 0.15)}

\item{batch_size}{Batch size (default: 64L)}
}
\value{
a list that include a training and a testing dataset, where both of them are 
				tfdataset object
}
\description{
Split a tfdataset object into training and testing sets
}
