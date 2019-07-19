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

//' Squared l2 Distance between two windows
//'
//' @description
//' This function gives distance matrix between two windows in multivariate time series
//' @param x NumericMatrix. first window.
//' @param y NumericMatrix. second window.
//' @return double. squared l2 distance matrix form.
//' @details
//' Compute
//' \deqn{d = x_{ij} - y_{kl}}
//' element-wise. After that,
//' \deqn{\frac{\sqrt{\sum{d}}}{wp}}
//' where w is the wize of window, and p is the number of variables.
//' @references
//' Yun, J.-H., Hwang, Y., Lee, W., Ahn, H.-K., & Kim, S.-K. (2018). \emph{Statistical Similarity of Critical Infrastructure Network Traffic Based on Nearest Neighbor Distances} (Vol. 11050, pp. 1–23). Presented at the Research in Attacks, Intrusions, and Defenses, Cham: Springer International Publishing. \url{http://doi.org/10.1007/978-3-030-00470-5_27}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double compute_euc(NumericMatrix x, NumericMatrix y) {
  int win = x.nrow();
  int px = x.ncol();
  double dist = 0;

  if (win != y.nrow())
    stop("x and y should have same column number");

  if (px != y.ncol())
    stop("x and y should have same row number");

  NumericVector l2(win * px);

  for (int i = 0; i < win; i++) {
    for (int j = 0; j < px; j++) {
      l2[i * win + j] = pow(x(i, j) - y(i, j), 2);
    }
  }

  dist = sqrt(sum(l2)) / (win * px);

  return dist;
}

// [[Rcpp::export]]
NumericMatrix sub_mat(NumericMatrix x, NumericVector row, NumericVector col) {
  NumericMatrix y(row.size(), col.size());

  for (int r = 0; r < row.size(); r++) {
    for (int c = 0; c < col.size(); c++) {
      y(r, c) = x(row[r], col[c]);
    }
  }

  return y;
}

// [[Rcpp::export]]
NumericVector partnnd(
  NumericMatrix data,
  int win,
  int jump,
  NumericVector partition,
  int base_id,
  Function d
) {
  NumericVector x_id = partition[partition == base_id];
  NumericVector y_id = partition[partition != base_id];

  int x_win = (x_id.size() - win) / jump + 1;
  NumericVector x_nnd(x_win);

  int y_win = (y_id.size() - win) / jump + 1;
  NumericVector y_nnd(y_win);

  NumericMatrix x(x_id.size(), data.ncol());
  NumericMatrix y(y_id.size(), data.ncol());

  for (int i = 0; i < x_win; i++) {
    for (int j = 0; j < y_win; j++) {
      y_nnd[j] = as<double>(d(x(Range(i * jump, i * jump + win - 1), _), y(Range(j * jump, j * jump + win - 1), _)));
    }
    x_nnd[i] = min(y_nnd);
  }

  return x_nnd;
}

//' Compute NND of normal data given window distance function
//'
//' @description
//' This function computes NND of a normal data set corresponding to defined window distance.
//' @param data NumericMatrix. data to be calculated NND.
//' @param part int. the number of partition.
//' @param win int. window size.
//' @param jump int. shift size.
//' @param d Function. distance function. By default, \code{\link{compute_euc}}. Arguments must be two matrices and the output double.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector. NND for each window. Its size is affected by "jump".
//' @details
//' First partition the whole given data set by "part".
//' Without the loss of generality, consider the first window.
//' We do not compute the distances between windows in the same window, assuming they are similar to each other.
//' Then for the first window, we consider every other slided window in the other partion.
//' Compute their distances. The minimum among the distance values is NND of the first window.
//' Repeat this procedure for every window in the partition.
//' Same for the other partition.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector nnd_normal(
  NumericMatrix data,
  int part,
  int win,
  int jump,
  Function d = compute_euc,
  bool display_progress = false
) {
  int n = data.nrow();
  int part_num = n / part;
  int part_rem = n % part;

  NumericMatrix x(win, data.ncol()); // window
  NumericMatrix y(win, data.ncol()); // versus window

  Progress p(part, display_progress);

  IntegerVector id = seq_len(part) - 1;
  IntegerVector idx = rep_each(id, part_num); // partition index
  IntegerVector idy(idx.size() + part_rem);
  idy[Range(0, idx.size() - 1)] = idx;
  if (part_rem != 0)
    idy[Range(idx.size(), idx.size() + part_rem - 1)] = rep(part - 1, part_rem);

  NumericVector nnd(idy.size());

  for (int i = 0; i < part; i++)
    nnd[Range(i * jump, i * jump + win - 1)] = partnnd(data, win, jump, idy, i, d);

  return nnd;
}

//' Compute NND for Online Data-set
//'
//' @description
//' This function computes NND of the data set which consists of anomaly.
//' @param data NumericMatrix. original data set.
//' @param newdata NumericMatrix. updated data set.
//' @param win int. window size.
//' @param jump int. shift size.
//' @param d Function. distance function. By default, \code{\link{compute_euc}}. Arguments must be two matrices and the output double.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector. NND for each window.
//' @details.
//' Measure a distance between windows in new data set and normal data set.
//' For this, slide windows in both set. Compute the distances and record them.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector pred_nnd(
  NumericMatrix data,
  NumericMatrix newdata,
  int win,
  int jump,
  Function d = compute_euc,
  bool display_progress = false
) {
  int n = data.nrow();
  NumericMatrix dat_new(n + win, data.ncol());
  NumericVector nnd_new((newdata.nrow() - win) / jump + 1);

  IntegerVector id(n + win);
  id[Range(0, n - 1)] = rep(0, n);
  id[Range(n, n + win - 1)] = rep(1, win);
}

