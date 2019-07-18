#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector rep_bool(bool x, int n) {
  LogicalVector y(n);

  for (int i = 0; i < n; i++) {
    y[i] = x;
  }

  return y;
}

// [[Rcpp::export]]
LogicalVector detect(NumericVector y, int win, int jump, double thr) {
  int n = y.size();
  int win_num = (n - win) / jump + 1;

  LogicalVector win_out(win);
  LogicalVector x(n);

  LogicalVector t_seq = rep_bool(true, win);
  LogicalVector f_seq = rep_bool(false, win);

  for (int i = 0; i < win_num; i++) {
    for (int j = 0; j < win; j++) {
      win_out[j] = y[i * jump + j] > thr;
    }

    if ( is_true(any(win_out)) ) {
      x[Range(i * jump, i * jump + win - 1)] = t_seq;
    } else {
      x[Range(i * jump, i * jump + win - 1)] = f_seq;
    }
  }

  return x;
}

//' Expand Label into Window
//'
//' @description
//' When detection label is given only one value per window,
//' this function expands it.
//'
//' Warning: Do not use jump parameter against confusion.
//' @param x LogicalVector detected label
//' @param win int Window size
//' @return LogicalVector
//' @details
//' In the window, the values are the same.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector expand_label(LogicalVector x, int win) {
  int win_num = x.size();
  int n = win_num * win;

  LogicalVector y(n);

  for (int i = 0; i < win_num; i++) {
    for (int j = 0; j < win; j++) {
      y[i * win + j] = x[i];
    }
  }

  return y;
}
