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

#' Anomaly detection after conducting NND
#'
#' @description
#' This function detects anomaly based on NND, given \code{\link{nns_cpp}}.
#' @param nnd NumericVector result of \code{\link{nns_cpp}}
#' @param win int window size for sliding window
#' @param thr threshold for anomaly detection, in each window
#' @param anom vector represents anomaly and normal. By default, TRUE and FALSE.
#' @return LogicalVector,
#' If NND is (strictly) larger than threshold then TRUE.
#' Otherwise, FALSE
#' @details
#' Given n x p data, slide a window.
#' Compute NND for each pair of moving window.
#' For threshold, users can use tail value of \code{\link{euc_pdf}}.
#' @references
#' Filonov, P., Kitashov, F., & Lavrentyev, A. (2017). \emph{RNN-based Early Cyber-Attack Detection for the Tennessee Eastman Process}. CoRR.
#'
#' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
#' @importFrom dplyr case_when
#' @export
detect_nndvec <- function(nnd, win, thr, anom = c(TRUE, FALSE)) {
  x <- nnd > thr
  result <-
    dplyr::case_when(
      x == TRUE ~ anom[1],
      x == FALSE ~ anom[2]
    )
  rep(result, each = win)
}
# LogicalVector detect_nndvec(NumericVector nnd, int win, double thr) {
#   int w = nnd.length();
#   LogicalVector x(w);
#
#   for (int i = 0; i < w; i++) {
#     x[i] = nnd[i] > thr;
#   }
#
#   return x;
# }
