#include <Rcpp.h>
using namespace Rcpp;

//' Matrix Representation of Sliding Windows
//'
//' @description Bind every window by row.
//' @param x NumericVector univariate data.
//' @param win int window size.
//' @param jump int jump size of sliding window.
//' @return NumericMatrix
//' @details
//' Each column represents an observation in each window.
//' Each row represents an distinct window.
//' Given the size of window win and the jump size jump, the number of the row becomes (nrow - win) / jump + 1.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericMatrix win_mat(NumericVector x, int win, int jump) {
  int win_num = (x.size() - win) / jump + 1;
  NumericMatrix out(win_num, win);

  for (int i = 0; i < win_num; i++)
    out(i, _) = x[Range(i * jump, i * jump + win - 1)];

  return out;
}
