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
//' @references Pavel Filonov, Andrey Lavrentyev, and Artem Vorontsov. 2016. \emph{Multivariate industrial time series with cyber-attack simulation: Fault detection using an lstm-based predictive data model}. arXiv preprint \url{arXiv:1612.06676} (2016).
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
