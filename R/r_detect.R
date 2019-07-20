#' Anomaly detection after conducting specific algorithm
#'
#' @description
#' This function detects anomaly based on result vector.
#' @param y NumericVector result of anomaly detection algorithm.
#' @param win int window size for sliding window.
#' @param jump int jump size for sliding window.
#' @param thr threshold for anomaly detection, in each window
#' @param label vector represents anomaly and normal. By default, TRUE and FALSE, respectively.
#' @return LogicalVector,
#' If NND is (strictly) larger than threshold then TRUE.
#' Otherwise, FALSE
#' @details
#' Large NND represents distinct pattern from the other windows.
#' Given pre-specified threshold \eqn{t}, we find the points such that
#' \deqn{d_i > t}
#' If at least one point in the window is determined to be an anomaly, this function will output the whole window as anomaly.
#' @seealso
#'     \code{\link{pred_nnd}}
#'     \code{\link{kl_fix}}
#'     \code{\link{kl_dynamic}}
#' @references
#' Filonov, P., Kitashov, F., & Lavrentyev, A. (2017). \emph{RNN-based Early Cyber-Attack Detection for the Tennessee Eastman Process}. CoRR.
#'
#' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
#' @importFrom dplyr case_when
#' @export
detect_anomaly <- function(y, win, jump, thr, label = c(TRUE, FALSE)) {
  x <- detect(y, win, jump, thr) # true and false
  dplyr::case_when(
    x == TRUE ~ label[1],
    x == FALSE ~ label[2]
  )
}
