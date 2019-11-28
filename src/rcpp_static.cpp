#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;
#include "misc.h"

//' Simple Static Threshold
//'
//' @description
//' This function detects anomaly for each window using static threshold.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param threshold double threshold for anomaly
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return LogicalVector
//' @details
//' If at least one observation in the window is larger than the threshold, the entire window is anomaly.
//'
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_static(NumericMatrix x, int win, int jump, double threshold, bool display_progress = false) {
  int n = x.nrow();
  int np = x.ncol();
  int win_num = (n - win) / jump + 1;

  LogicalVector anomaly(win_num);
  NumericMatrix batch(win, np);

  Progress p(win_num - 1, display_progress);
  for (int i = 0; i < win_num; i++) {

    if (Progress::check_abort())
      return -1.0;
    p.increment();

    batch = x(Range(i * jump, i * jump + win - 1), _);

    anomaly[i] = all(batch > threshold).is_true();
  }

  return anomaly;
}

// [[Rcpp::export]]
NumericVector compute_mse(NumericMatrix x, int win, int jump) {
  int n = x.nrow();
  int np = x.ncol();
  int win_num = (n - win) / jump + 1;

  NumericMatrix batch(win, np);
  NumericVector mse(win_num);
  NumericVector mse_win(win);

  for (int i = 0; i < win_num; i++) {
    batch = x(Range(i * jump, i * jump + win - 1), _);
    for (int j = 0; j < win; j++) {
      mse_win[j] = sum(pow(batch(j, _), 2.0)) / np;
    }
    mse[i] = min(mse_win);
  }

  return mse;
}

//' Static Threshold based on MSE
//'
//' @description
//' This function detects anomaly for each window based on window MSE.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param threshold double threshold for anomaly
//' @return LogicalVector
//' @details
//' If the MSE is larger than given threshold, then the window is anomaly.
//'
//' @references Pavel Filonov, Andrey Lavrentyev, and Artem Vorontsov. 2016. \emph{Multivariate industrial time series with cyber-attack simulation: Fault detection using an lstm-based predictive data model}. arXiv preprint \url{arXiv:1612.06676} (2016).
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_mse(NumericMatrix x, int win, int jump, double threshold) {
  int n = x.nrow();
  int win_num = (n - win) / jump + 1;

  LogicalVector anomaly(win_num);
  NumericVector mse(win_num);
  mse = compute_mse(x, win, jump);

  anomaly = mse > threshold;

  return anomaly;
}

//' p-norm of each observation
//'
//' @description
//' This function computes p-norm w.r.t. the variables.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param norm int p-norm
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector
//' @seealso
//'   \code{\link{detect_norm}}
//'   \code{\link{train_norm}}
//'
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector compute_norm(NumericMatrix x, int norm, bool display_progress = false) {
  int n = x.nrow();
  NumericVector error(n);
  // p-norm
  Progress p(n - 1, display_progress);
  for (int i = 0; i < n; i++) {
    if (Progress::check_abort())
      return -1.0;
    p.increment();

    error[i] = sum(pow(x(i, _), norm));
  }
  return error;
}

//' Static Threshold based on p-norm
//'
//' @description
//' This function detects anomaly for each observation based on p-norm.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param norm int p-norm
//' @param threshold double threshold for anomaly
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return LogicalVector
//' @details
//' If the p-norm is larger than given threshold, then the observation is anomaly.
//'
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_norm(NumericMatrix x, int norm, double threshold, bool display_progress = false) {
  int n = x.nrow();

  LogicalVector anomaly(n);
  NumericVector error(n);
  error = compute_norm(x, norm, display_progress);
  anomaly = error > threshold;

  return anomaly;
}

//' Window CUSUM
//'
//' @description
//' This function computes the sum of p-norm in each window.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param norm int p-norm
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector
//' @details
//' This function will be used in both \code{\link{detect_cusum}} and \code{\link{train_cusum}}
//' Use error data as an input.
//'
//' @seealso
//'   \code{\link{detect_cusum}}
//'   \code{\link{train_cusum}}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector compute_cusum(NumericMatrix x, int win, int jump, int norm, bool display_progress = false) {
  int n = x.nrow();
  int win_num = (n - win) / jump + 1;
  LogicalVector anomaly(win_num);
  NumericVector error(win);
  NumericVector error_sum(win_num);

  Progress p(win_num - 1, display_progress);
  for (int i = 0; i < win_num; i++) {
    if (Progress::check_abort())
      return -1.0;
    p.increment();

    for (int w = 0; w < win; w++) {
      error[w] = sum(pow(x(i * jump + w, _), norm));
    }

    error_sum[i] = sum(error);
  }
  return error_sum;
}

//' CUSUM
//'
//' @description
//' This function detects anomaly for each window based on CUSUM.
//'
//' @param x NumericMatrix multivariate time series, which is forecasting error
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param norm int p-norm
//' @param threshold double threshold for anomaly
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return LogicalVector
//' @details
//' If the p-norm is larger than given threshold, then the observation is anomaly.
//'
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_cusum(NumericMatrix x, int win, int jump, int norm, double threshold, bool display_progress = false) {
  int n = x.nrow();
  int win_num = (n - win) / jump + 1;
  LogicalVector anomaly(win_num);
  NumericVector error_sum(win_num);
  error_sum = compute_cusum(x, win, jump, norm, display_progress);
  anomaly = error_sum > threshold;
  return anomaly;
}

