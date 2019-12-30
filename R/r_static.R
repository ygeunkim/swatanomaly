#' Choose Threshold for p-norm
#'
#' @param x Train data
#' @param norm p-norm
#' @param prob Quantile
#' @return Suggested threshold value for \code{\link{detect_norm}}
#' @details
#' By computing the quantile of every p-norm, threshold is given.
#' @seealso
#' \code{\link{compute_norm}}
#' \code{\link{detect_norm}}
#' @importFrom stats quantile
#' @export
train_norm <- function(x, norm, prob = .9) {
  if (!is.matrix(x)) x <- as.matrix(x)
  error <- compute_norm(x, norm)
  quantile(error, probs = prob)
}

#' Choose Threshold for CUSUM
#'
#' @param x Train data
#' @param win window size
#' @param jump jump size
#' @param norm p-norm
#' @param prob Quantile
#' @return Suggested threshold value for \code{\link{detect_cusum}}
#' @details
#' By computing the quantile of every CUSUM, threshold is given.
#' @seealso
#' \code{\link{compute_cusum}}
#' \code{\link{detect_cusum}}
#' @importFrom stats quantile
#' @export
train_cusum <- function(x, win, jump, norm, prob = .9) {
  if (!is.matrix(x)) x <- as.matrix(x)
  error <- compute_cusum(x, win, jump, norm)
  quantile(error, probs = prob)
}
