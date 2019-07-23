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
NumericMatrix sub_mat(NumericMatrix x, IntegerVector row, IntegerVector col) {
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

// [[Rcpp::export]]
NumericVector concat_vec(NumericVector x, NumericVector y) {
  int nx = x.size();
  int ny = y.size();
  NumericVector out(nx + ny);

  out[Range(0, nx - 1)] = x;
  out[Range(nx, ny - 1)] = y;

  return out;
}

// [[Rcpp::export]]
double compute_q7(NumericVector x, double prob) {
  NumericVector y = sort_unique(x);
  int n = y.size();

  int h = n * (prob - 1e-9);

  return y[h];
}

// [[Rcpp::export]]
double nrd0(NumericVector x) {
  if (x.size() < 2)
    stop("need at least 2 data points");

  double hi = sd(x);
  double lo = std::min(hi, (compute_q7(x, .75) - compute_q7(x, .25)) / 1.34);

  if (lo == 0)
    (lo = hi) || (lo = abs(x[0])) || (lo = 1);

  return .9 * lo * pow(x.size(), -.2);
}



