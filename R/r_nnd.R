#' NND pdf Threshold
#'
#' @description
#' This function computes threshold using \code{\link{nnd_normal}}
#' @param data NumericMatrix multivariate data set
#' @param part_size part for \code{\link{nnd_normal}}
#' @param win_size win for \code{\link{nnd_normal}}
#' @param jump_size jump for \code{\link{nnd_normal}}
#' @param prob probability for quantile. See \link[stats]{quantile.default}.
#' @param display_progress If TRUE, display a progress bar. By default, FALSE.
#' @return double quantile value for the result of \code{\link{nnd_normal}}
#' @importFrom stats quantile
#' @details
#' Sometimes you might use NND pdf for threshold.
#' Threshold is a tail of pdf, e.g. 0.99.
#' @seealso
#'    \code{\link{nnd_normal}}
#'    \link[stats]{quantile.default}
#'    \code{\link{detect_anomaly}}
#' @export
nnd_thr <- function(data, part_size, win_size, jump_size, prob, display_progress = FALSE) {
  quantile(
    nnd_normal(data, part_size, win_size, jump_size, display_progress),
    probs = prob
  )
}
