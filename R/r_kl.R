#' Gaussian Kernel Density Estimation
#'
#' @description
#' This function calls \link[stats]{density.default}.
#' @param x NumericVector data for estimation. x of \link[stats]{density.default}.
#' @param ... Additional arguments of \link[stats]{density.default}.
#' @return NumericMatrix of 2 columns.
#' First column is the n coordinates of the points where the density is estimated (x).
#' Second column is the estimated density values (y).
#' @details
#' Note that \link[stats]{density.default} has various arguments and options.
#' This function, howerver, only computes gaussian kernel.
#' It tries to detect windows derived from other Normal distribution.
#' @seealso
#'    \link[stats]{density.default}
#'    \code{\link{density_cpp}}
#' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
#' @importFrom stats density
#' @export
est_density <- function(x, ...) {
  den <- density(x, ...)
  matrix(
    c(den$x, den$y),
    ncol = 2,
    dimnames = list(NULL, c("x", "y"))
  )
}
