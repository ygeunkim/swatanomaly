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
LogicalVector rep_bool(bool x, int n) {
  LogicalVector y(n);

  for (int i = 0; i < n; i++) {
    y[i] = x;
  }

  return y;
}

// [[Rcpp::export]]
NumericMatrix rbind_mat(NumericMatrix x, NumericMatrix y) {
  int nx = x.nrow();
  int ny = y.nrow();
  NumericMatrix out(nx + ny, x.ncol());

  for (int i = 0; i < nx; i++)
    out(i, _) = x(i, _);

  for (int j = 0; j < ny; j++)
    out(nx + j, _) = y(j, _);

  return out;
}




