% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_reads.r
\name{count_reads,GRanges,character,BSgenome-method}
\alias{count_reads,GRanges,character,BSgenome-method}
\title{Count reads}
\usage{
\S4method{count_reads}{GRanges,character,BSgenome}(x, filenames, genome, fragment_size_range = c(0L, 320))
}
\arguments{
\item{x}{a GRanges object defining the peaks; the width must be the same for all peaks.}

\item{filenames}{BAM file name(s).}

\item{genome}{a BS genome object such as BSgenome.Mmusculus.UCSC.mm10}

\item{fragment_size_range}{fragment size ranges (default:  c(80, 320))}
}
\value{
read counts
}
\description{
Count how many reads fall into a specific fragment size range
}
\author{
Wuming Gong (gongx030@umn.edu)
}
