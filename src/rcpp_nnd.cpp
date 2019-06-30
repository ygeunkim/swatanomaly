#include <Rcpp.h>
using namespace Rcpp;

//' Sums of squares in C++
//'
//' @description Compute a SS in C++
//' @param x NumericVector
//' @return double
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double sum_sq(NumericVector x) {
  int n = x.size();
  double sum = 0;

  for (int i = 0; i < n; i ++) {
    sum += pow(x[i], 2.0);
  }

  return sum;
}

//' Euclidean pdf in Rcpp
//'
//' @description
//' This function computes a euclidean NND pdf of two multivariate series using Rcpp. See details what it is.
//' @param x NumericMatrix column should indicate variable
//' @param y NumericMatrix column should indicate variable
//' @param partition int partitioning equally the series
//' @return NumericVector NND for each block
//' @details
//' First partitioning the series equally.
//' Next for each partitioned block, it calculates sqrt(sum((x_i - y_i)^2)) versus the other blocks.
//' Find the minimum result for each block. This is NND of each block.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector euc_pdf(NumericMatrix x, NumericMatrix y, int partition) {
  NumericVector nnd(partition);
  NumericVector euc(partition);
  int nx = x.nrow();
  int ny = y.nrow();

  for (int i = 0; i < nx; i += (int)partition) {
    for (int j = 0; j < ny; j += (int)partition) {
      for (int k = 0; k < partition; k ++) {
        nnd[k] = sqrt(sum_sq(x(i, i + partition - 1) - y(j, j + partition - 1)));
      }
      euc[i] = min(nnd);
    }
  }

  return euc;
}

//' Euclidean distance between two matrices in Rcpp
//'
//' @description
//' This function computes a euclidean distance between two multivariate series using Rcpp.
//' @param x NumericMatrix column should indicate variable
//' @param y NumericMatrix column should indicate variable
//' @return double
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double euc_dist(NumericMatrix x, NumericMatrix y) {
  int nx = x.nrow();
  int px = x.ncol();
  int ny = y.nrow();
  int py = y.ncol();
  double euc = 0;

  if (nx != ny | px != py) {
    stop("x and y should be have same dimension");
  }

  for (int i = 0; i < nx; i++) {
    euc += sum_sq(x(i, _) - y(i, _));
  }

  euc = sqrt(euc);
  return euc;
}

//' Windowed NNS in Rcpp
//'
//' @description
//' This function computes a windowed NNS.
//' Compute NND sliding window across given series.
//' @param data NumericMatrix multivariate data set
//' @param win int window size for sliding window
//' @return NumericMatrix
//' @details
//' Given n x p data, slide a window.
//' Compute NND for each pair of moving window.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector nns_cpp(NumericMatrix data, int win) {
  int n = data.nrow();
  NumericVector x1(win);
  NumericVector x2(win);
  NumericMatrix distmat(n + win - 1, n + win - 1);
  NumericVector distvec(n + win - 1);
  NumericVector findmin(n + win - 1);

  for (int i = 0; i < n + win - 1; i++) {
    for (int j = 0; j < n + win - 1; j++) {
      distmat(i, j) = euc_dist(data(Range(i, i + win - 1), _), data(Range(j, j + win - 1), _));
    }
    findmin = distmat(i, _);
    findmin[i] = max(findmin);
    distvec[i] = min(findmin);
  }

  return distvec;
}

//' Anomaly detection using NND
//'
//' @description
//' This function detects anomaly based on NND.
//' @param data NumericMatrix multivariate data set
//' @param win int window size for sliding window
//' @param thr threshold for anomaly detection, in each window
//' @return LogicalVector,
//' If NND is (strictly) larger than threshold then TRUE.
//' Otherwise, FALSE
//' @details
//' Given n x p data, slide a window.
//' Compute NND for each pair of moving window.
//' For threshold, users can use tail value of \code{\link{euc_pdf}}.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_nnd(NumericMatrix data, int win, NumericVector thr) {
  NumericVector distvec = nns_cpp(data, win);
  int w = thr.length();
  LogicalVector x(w);

  for (int i = 0; i < w; i++) {
    x[i] = distvec[i] > thr[i];
  }

  return x;
}
