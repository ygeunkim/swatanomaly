#' NND pdf Threshold
#'
#' @description
#' This function computes threshold using \code{\link{euc_pdf}}
#' @param data NumericMatrix multivariate data set
#' @param part_size partition for \code{\link{euc_pdf}}
#' @param prob probability for quantile. See \link[stats]{quantile.default}.
#' @param display_progress If TRUE, display a progress bar. By default, FALSE.
#' @return double quantile value for the result of \code{\link{euc_pdf}}
#' @importFrom stats quantile
#' @details
#' Sometimes you might use NND pdf for threshold.
#' Threshold is a tail of pdf, e.g. 0.99.
#' @seealso
#'    \code{\link{euc_pdf}}
#'    \link[stats]{quantile.default}
#'    \code{\link{detect_nnd}}
#'    \code{\link{detect_nndvec}}
#' @export
nnd_thr <- function(data, part_size, prob, display_progress) {
  quantile(
    euc_pdf(data, part_size, display_progress),
    probs = prob
  )
}
