% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{detect_norm}
\alias{detect_norm}
\title{Static Threshold based on p-norm}
\usage{
detect_norm(x, norm, threshold)
}
\arguments{
\item{x}{NumericMatrix multivariate time series, which is forecasting error}

\item{norm}{int p-norm}

\item{threshold}{double threshold for anomaly}
}
\value{
LogicalVector
}
\description{
This function detects anomaly for each observation based on p-norm.
}
\details{
If the p-norm is larger than given threshold, then the observation is anomaly.
}
