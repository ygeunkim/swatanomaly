// #include <Rcpp.h>
#include <armadillo>
#include <RcppArmadillo.h>
// using namespace Rcpp;
using namespace arma;

//' Matrix Representation of Sliding Windows using PCA
//'
//' @description
//' Bind every window by row,
//' using principal components analysis in each window.
//' @param x mat multivariate series.
//' @param win int window size.
//' @param jump int jump size of sliding window.
//' @return mat
//' @details
//' Each column represents an observation in each window.
//' Each row represents an distinct window.
//' Given the size of window win and the jump size jump, the number of the row becomes (nrow - win) / jump + 1.
//' The data set is multivariate series, in general, so we perform dimension reduction using PCA.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::mat reduce_pca(arma::mat x, int win, int jump) {
  int win_num = (x.n_rows - win) / jump + 1;

  arma::mat out(win_num, win);
  arma::mat z(win, x.n_cols);

  for (int i = 0; i < win_num; i++) {
    z = princomp(x.rows(i * jump, i * jump + win - 1));
    out.row(i) = z.col(0);
  }

  return out;
}

//' Matrix Representation of Sliding Windows using MDS
//'
//' @description
//' Bind every window by row,
//' using multidimensional scaling in each window.
//' @param x mat multivariate series.
//' @param win int window size.
//' @param jump int jump size of sliding window.
//' @return mat
//' @details
//' Each column represents an observation in each window.
//' Each row represents an distinct window.
//' Given the size of window win and the jump size jump, the number of the row becomes (nrow - win) / jump + 1.
//' The data set is multivariate series, in general, so we perform dimension reduction using MDS.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
arma::mat reduce_mds(arma::mat x, int win, int jump) {
  int win_num = (x.n_rows - win) / jump + 1;

  arma::mat out(win_num, win);
  arma::mat z(win, x.n_cols);
  arma::mat b = x * x.t(); // subset i * jump x i * jump to i * jump + win - 1 x i * jump + win - 1 in the loop

  arma::mat l(win, win);

  arma::vec eigenval;
  arma::mat eigenvec;

  for (int i = 0; i < win_num; i++) {
    eig_sym(eigenval, eigenvec, x.submat(i * jump, i * jump, i * jump + win - 1, i * jump + win - 1));
    l = sqrt(diagmat(eigenval));
    out.row(i) = eigenvec.col(0) * l.col(0);
  }

  return out;
}
