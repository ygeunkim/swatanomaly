% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{partnnd}
\alias{partnnd}
\title{Sliding window for NND}
\usage{
partnnd(data, win, jump, from, to)
}
\arguments{
\item{data}{NumericMatrix. data to be calculated NND.}

\item{win}{int. window size.}

\item{jump}{int. shift size.}

\item{from}{int. First index of chosen partition.}

\item{to}{int. Last index of chosen partition.}
}
\value{
NumericVector. NND vector for each window in the chosen partition.
}
\description{
For the chosen partition, get NND for every pair of window.
This function is constructed for the other function.
}
\seealso{
\code{\link{nnd_normal}}
    \code{\link{pred_nnd}}
}
