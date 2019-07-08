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

//' Euclidean NND between validation and training blocks
//'
//' @description
//' This function computes NND corresponding to euclidean distance between validation sets and training sets built by cross-validation.
//' @param x NumericMatrix validation set.
//' @param y NumericMatrix training set.
//' @return NumericVector length identical to the validation set
//' @details
//' Consider input x and y. Compute
//' \deqn{\sqrt{\sum (x_{ij} - y_{ij})^2}}
//' for each observation.
//' Calculate Euclidean distance between a point in validation series versus each point in training series.
//' Find minimum value. It is NND.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector euc_nnd(NumericMatrix x, NumericMatrix y) {
  int nx = x.nrow();
  int px = x.ncol();
  int ny = y.nrow();
  int py = y.ncol();

  NumericVector euc(ny);
  NumericVector nnd(nx);

  if (px != py) {
    stop("x and y should be have same dimension");
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      euc[j] = sum_sq(x(i, _) - y(j, _));
    }
    nnd[i] = min(euc);
  }

  nnd = sqrt(nnd);
  return nnd;
}

//' Remove row index of a matrix in C++
//'
//' @description
//' This function removes a row index of NumericMatrix in Rcpp.
//' @param x NumericMatrix
//' @param rowID IntegerVector row ids to be removed.
//' @return NumericMatrix
//' @useDynLib swatanomaly
//' @references \url{https://stackoverflow.com/questions/33507695/rcpp-numericmatrix-how-to-erase-a-row-column}
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericMatrix row_erase(NumericMatrix x, IntegerVector rowID) {
  // rowID = rowID.sort();

  NumericMatrix x2(Dimension(x.nrow() - rowID.size(), x.ncol()));

  int iter = 0;
  int del = 1; // to count deleted elements

  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}

//' Sequence by 1 in Rcpp
//'
//' @description
//' This function generates a integer sequence with increment of 1 in Rcpp.
//' @param from int the starting value of the sequence.
//' @param to int the end value of the sequence.
//' @return IntegerVector
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
IntegerVector seq_rcpp(int from, int to) {
  IntegerVector x(to - from + 1);

  for (int i = 0; i < x.size(); i++) {
    x[i] = from + i;
  }

  return x;
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
//' Implement k-fold cross-validation. Here, k is partition.
//' First partitioning the series equally.
//' One block is validation set, and the others are trainig.
//' Next for each validation series, it calculates sqrt(sum((x_i - x_j)^2)) versus training series.
//' Find the minimum result for each validation block. This is NND of each block.
//' Finally, you can get NND for every block and this is pdf for NND.
//' For \code{\link{detect_nndvec}}, this pdf is able to threshold.
//' Threshold is a tail of pdf, e.g. 0.99.
//' @seealso
//'    \code{\link{euc_dist}}
//'    \code{\link{nnd_thr}}
//'    \code{\link{detect_nndvec}}
//' @references
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector euc_pdf(NumericMatrix x, int partition, bool display_progress = false) {
  int div = x.nrow() / partition;

  NumericVector euc(x.nrow());

  Progress p(partition, display_progress);

  for (int i = 0; i < partition; i++) {

    if (Progress::check_abort())
      return -1.0;

    p.increment();

    euc[Range(i * div, i * div + div - 1)] = euc_nnd(x(Range(i * div, i * div + div - 1), _), row_erase(x, seq_rcpp(i * div, i * div + div - 1)));
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
//' Given n x p data, partition a window.
//' Compute NND for each pair of window.
//' The method is similar to \code{\link{euc_pdf}}.
//' Note that the number of windows is nrow - win + 1 given size of window, win.
//' @references
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector nns_cpp(NumericMatrix data, int win, bool display_progress = false) {
  int n = data.nrow();
  int win_num = n / win;

  NumericVector sliding(win_num);
  NumericVector distvec(n);

  Progress p(win_num, display_progress);

  for (int i = 0; i < win_num; i++) {

    if (Progress::check_abort())
      return -1.0;

    p.increment();

    distvec[Range(i * win, i * win + win - 1)] = euc_nnd(data(Range(i * win, i * win + win - 1), _), row_erase(data, seq_rcpp(i * win, i * win + win - 1)));
  }

  return distvec;
}

// [[Rcpp::export]]
LogicalVector rep_bool(bool x, int n) {
  LogicalVector y(n);

  for (int i = 0; i < n; i++) {
    y[i] = x;
  }

  return y;
}

// [[Rcpp::export]]
LogicalVector detect_nnd(NumericVector nnd, int win, double thr) {
  int n = nnd.size();
  int win_num = n / win;

  LogicalVector win_out(win);
  LogicalVector x(n);

  LogicalVector t_seq = rep_bool(true, win);
  LogicalVector f_seq = rep_bool(false, win);

  for (int i = 0; i < win_num; i++) {
    for (int j = 0; j < win; j++) {
      win_out[j] = nnd[i * j + j] > thr;
    }

    if ( is_true(any(win_out)) ) {
      x[Range(i * win, i * win + win - 1)] = t_seq;
    } else {
      x[Range(i * win, i * win + win - 1)] = f_seq;
    }
  }

  return x;
}

