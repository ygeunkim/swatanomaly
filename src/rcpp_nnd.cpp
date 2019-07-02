#include <Rcpp.h>
#include <progress.hpp>
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

//' Euclidean distance between two matrices in Rcpp
//'
//' @description
//' This function computes a euclidean distance between two multivariate series using Rcpp.
//' @param x NumericMatrix column should indicate variable
//' @param y NumericMatrix column should indicate variable
//' @return double
//' @details
//' For input x and y, compute
//' \deqn{\sum \sqrt{\sum (x_{ij} - y_{ij})^2}}
//' At first, the function calculates Euclidean distance pairwisely.
//' After that, sum over every observation.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double euc_dist(NumericMatrix x, NumericMatrix y) {
  int nx = x.nrow();
  int px = x.ncol();
  int ny = y.nrow();
  int py = y.ncol();
  NumericVector euc(px);

  if (nx != ny | px != py) {
    stop("x and y should be have same dimension");
  }

  for (int j = 0; j < px; j++) {
    euc[j] = sum_sq(x(_, j) - y(_, j));
  }

  double dist = sum(sqrt(euc));
  return dist;
}

//' Euclidean pdf in Rcpp
//'
//' @description
//' This function computes a euclidean NND pdf of two multivariate series using Rcpp. See details what it is.
//' @param x NumericMatrix, column should indicate variable.
//' @param partition int, equally partitioning the series.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector, NND pdf
//' @details
//' First partitioning the series equally.
//' Next for each partitioned block, it calculates sqrt(sum((x_i - x_j)^2)) versus the other partioned ones.
//' Find the minimum result for each block. This is NND of each block.
//' Finally, you can get NND for every block and this is pdf for NND.
//' For \code{\link{detect_nnd}}, this pdf is able to threshold.
//' Threshold is a tail of pdf, e.g. 0.99.
//' @seealso
//'    \code{\link{euc_dist}}
//'    \code{\link{nnd_thr}}
//'    \code{\link{detect_nnd}}
//' @references
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector euc_pdf(NumericMatrix x, int partition, bool display_progress = false) {
  int nx = x.nrow();
  int div = nx / partition;

  NumericVector nnd(div);
  NumericVector euc(div);

  Progress p(div * div, display_progress);

  for (int i = 0; i < div; i++) {

    if (Progress::check_abort())
      return -1.0;

    for (int j = 0; j < div; j ++) {
      p.increment();
      nnd[j] = euc_dist(x(Range(i, i + partition - 1), _), x(Range(j, j + partition - 1), _));
    }
    nnd[i] = max(nnd);
    euc[i] = min(nnd);
  }

  return euc;
}

//' Windowed NNS in Rcpp
//'
//' @description
//' This function computes a windowed NNS.
//' Compute NND sliding window across given series.
//' @param data NumericMatrix multivariate data set
//' @param win int window size for sliding window
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector, NND for each window index (index represented by its starting point)
//' @details
//' Given n x p data, slide a window.
//' Compute NND for each pair of moving window.
//' Note that the number of windows is nrow - win + 1 given size of window win.
//' @references
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector nns_cpp(NumericMatrix data, int win, bool display_progress = false) {
  int n = data.nrow();
  NumericVector sliding(n - win + 1);
  NumericVector distvec(n - win + 1);

  Progress p((n - win + 1) * (n - win + 1), display_progress);

  for (int i = 0; i < n - win + 1; i++) {

    if (Progress::check_abort())
      return -1.0;

    for (int j = 0; j < n - win + 1; j++) {
      p.increment();
      sliding[j] = euc_dist(data(Range(i, i + win - 1), _), data(Range(j, j + win - 1), _));
    }
    sliding[i] = max(sliding);
    distvec[i] = min(sliding);
  }

  return distvec;
}

//' Anomaly detection using NND
//'
//' @description
//' This function detects anomaly based on NND.
//' @param data NumericMatrix multivariate data set.
//' @param win int window size for sliding window.
//' @param thr double threshold that will be compared to nnd vector.
//' @return LogicalVector,
//' If NND is (strictly) larger than threshold then TRUE.
//' Otherwise, FALSE
//' @details
//' Given n x p data, slide a window.
//' Compute NND for each pair of moving window.
//' For threshold, users can use tail value of \code{\link{euc_pdf}}.
//' @references
//' Filonov, P., Kitashov, F., & Lavrentyev, A. (2017). \emph{RNN-based Early Cyber-Attack Detection for the Tennessee Eastman Process}. CoRR.
//'
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_nnd(NumericMatrix data, int win, double thr) {
  int n = data.nrow();
  NumericVector distvec = nns_cpp(data, win);

  LogicalVector x(n - win + 1);

  for (int i = 0; i < n - win + 1; i++) {
    x[i] = distvec[i] > thr;
  }

  return x;
}

//' Anomaly detection after conducting NND
//'
//' @description
//' This function detects anomaly based on NND, given \code{\link{nns_cpp}}.
//' @param nnd NumericVector result of \code{\link{nns_cpp}}
//' @param win int window size for sliding window
//' @param thr threshold for anomaly detection, in each window
//' @return LogicalVector,
//' If NND is (strictly) larger than threshold then TRUE.
//' Otherwise, FALSE
//' @details
//' Given n x p data, slide a window.
//' Compute NND for each pair of moving window.
//' For threshold, users can use tail value of \code{\link{euc_pdf}}.
//' @references
//' Filonov, P., Kitashov, F., & Lavrentyev, A. (2017). \emph{RNN-based Early Cyber-Attack Detection for the Tennessee Eastman Process}. CoRR.
//'
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector detect_nndvec(NumericVector nnd, int win, double thr) {
  int w = nnd.length();
  LogicalVector x(w);

  for (int i = 0; i < w; i++) {
    x[i] = nnd[i] > thr;
  }

  return x;
}

