#include <Rcpp.h>
using namespace Rcpp;

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
