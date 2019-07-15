#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;

//' Aggregate Multivariate Time Series for K-L divergence
//'
//' @description
//' This functions aggregates multivariate times series into univariate time series.
//' See details.
//' @param x NumericMatrix multivariate time series
//' @return NumericVector
//' @details
//' To eliminate local spikes, compute usual distance.
//' \deqn{\sum_{i \neq j}^p \lvert x_{ti} - x_{tj}}
//' where p is the number of variables, and t is the index of time.
//' This enables to explain the correlation between the series.
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector aggregate_mts(NumericMatrix x) {
  int p = x.ncol();
  NumericVector y(x.nrow());

  for (int i = 0; i < p; i++) {
    for (int j = 0; j < i; j++) {
      y += abs(x(_, i) - x(_, j));
    }
  }

  return y;
}

//' Gaussian Kernel Density Estimation in C++
//'
//' @description
//' This function calls \link[stats]{density.default} in Rcpp.
//' Since it calls an R function to Rcpp, it might be slower even than R function itself.
//' This one is just for the sake of defining other functions in Rcpp syntax more easily.
//' @param x NumericVector data for estimation. x of \link[stats]{density.default}.
//' @return NumericMatrix of 2 columns.
//' First column is the n coordinates of the points where the density is estimated (x).
//' Second column is the estimated density values (y).
//' @details
//' Note that \link[stats]{density.default} has various arguments and options.
//' This function, howerver, only computes gaussian kernel.
//' It tries to detect windows derived from other Normal distribution.
//' @seealso
//'    \link[stats]{density.default}
//'    \code{\link{est_density}}
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @importFrom stats density
//' @export
// [[Rcpp::export]]
NumericMatrix density_cpp(NumericVector x) {
  Function dens("density.default");

  List den = dens(x);
  if (! den.inherits("density")) stop("not an appropriate S3 density class");

  NumericMatrix den_xy(512, 2); // default length 512

  NumericVector den_x = den["x"];
  NumericVector den_y = den["y"];

  den_xy(_, 0) = den_x;
  den_xy(_, 1) = den_y;

  return den_xy;
}

//' Kullback-Leibler divergence estimation between two densities
//'
//' @description
//' This function computes Kullback-Leibler divergence from f2 to f1.
//' When sliding windows, f1 is the previous pdf and f2 is the current pdf.
//' @param f1 NumericMatrix density estimated by \code{\link{est_density}}. previous pdf.
//' @param f2 NumericMatrix density estimated by \code{\link{est_density}}. current pdf.
//' @return double
//' @details
//' Let \eqn{\mathcal{X}} be the support of f1. Then K-L divergence from f2 to f1 is defined by
//' \deqn{E_{X_1} \log \frac{f_1 (x)}{f_2 (x)}}
//' Probability mass is estimated from density by
//' \deqn{f \Delta x}
//' In turn, we can compute K-L divergence using mass p and q by
//' \deqn{\sum_{\mathcal{X}} p(x) \log \frac{p (x)}{q (x)}}
//' To estimate this value, first use \code{\link{est_density}} or \code{\link{density_cpp}}, and estimate density of each window.
//' @seealso
//'     \code{\link{est_density}}
//'     \code{\link{density_cpp}}
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double compute_kl(NumericMatrix f1, NumericMatrix f2) {
  double kl = 0;
  double sum_kl = 0;

  double px = 0;
  double qx = 0;

  for (int i = 0; i < f2.nrow() - 1; i++) {
    if (is_true( all(f1(_, 0) < f2(i, 0)) )) {
      sum_kl = 0;
    } else {
      qx = f1(i, 1) * abs(f2(i + 1, 0) - f2(i, 0));
      px = f2(i, 1) * abs(f2(i + 1, 0) - f2(i, 0));
      if (px == 0 | qx == 0) {
        sum_kl = 0;
      } else {
        sum_kl = px * log(px / qx);
      }
    }
    kl += sum_kl;
  }

  return kl;
}

//' Fixed lambda Algorithm
//'
//' @description
//' This function implement dynamic \eqn{\lambda} algorithm for K-L divergence based on gaussianity.
//' @param x NumericVector univariate data set.
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param lambda double threshold of K-L divergence.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector
//' @details
//' Slide windows. In each window, estimate a density based on \link[stats]{density.default}.
//' Between two neighbored window, K-L divergence from current density to previous density can be computed using \code{\link{compute_kl}}.
//' If K-L divergence is less than \eqn{\lambda}, update the threshold by
//' @seealso
//'    \link[stats]{density.default}
//'     \code{\link{est_density}}
//'     \code{\link{density_cpp}}
//'     \code{\link{compute_kl}}
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector kl_fix(NumericVector x, int win, int jump, double lambda, bool display_progress = false) {
  int n = x.size();
  int win_num = (n - win) / jump + 1;
  double kl_s = 0;

  if (lambda <= 0) stop("lambda should be positive");

  NumericMatrix f1(512, 2);
  NumericMatrix f2(512, 2);
  NumericVector kl(win_num - 1);

  Progress p(win_num - 1, display_progress);

  for (int i = 0; i < win_num - 1; i++) {
    if (Progress::check_abort())
      return -1.0;

    p.increment();

    if (kl_s < lambda) {
      f1 = density_cpp(x[Range(i * jump, i * jump + win - 1)]);
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2);
      kl_s = kl[i];
    } else {
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2); // f1 = of normal
      kl_s = kl[i];
    }
  }

  return kl;
}


//' Dynamic lambda Algorithm
//'
//' @description
//' This function implement dynamic \eqn{\lambda} algorithm for K-L divergence based on gaussianity.
//' @param x NumericVector univariate data set.
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param lambda_p double initializing lambda_p for the threshold.
//' @param eps double initializing epsilon for the threshold.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return NumericVector
//' @details
//' Basically, this algorithm use neighboring-window method.
//' Slide windows. In each window, estimate a density based on \link[stats]{density.default}.
//' Between two neighbored window, K-L divergence from current density to previous density can be computed using \code{\link{compute_kl}}.
//' Given \eqn{\lambda^{\prime}} and \eqn{\epsilon}, set threshold \eqn{\lambda} by \eqn{\lambda = \lambda^{\prime} \epsilon}.
//' If K-L divergence is less than \eqn{\lambda}, update the threshold by
//' \deqn{\lambda = \lambda^{\prime} (d_{j - 2} + \epsilon)}
//' where \eqn{d_{j - 2}} is the K-L divergence two-step before.
//' Otherwise, keep using the old one.
//' @seealso
//'    \link[stats]{density.default}
//'     \code{\link{est_density}}
//'     \code{\link{density_cpp}}
//'     \code{\link{compute_kl}}
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector kl_dynamic(NumericVector x, int win, int jump, double lambda_p, double eps, bool display_progress = false) {
  int n = x.size();
  int win_num = (n - win) / jump + 1;
  NumericMatrix f1(512, 2);
  NumericMatrix f2(512, 2);

  double lambda = lambda_p * eps;
  double kl_s = 0;

  Progress p(win_num - 3, display_progress);

  NumericVector kl(win_num - 1);

  // i = 1
  f1 = density_cpp(x[Range(0, win - 1)]);
  f2 = density_cpp(x[Range(jump, jump + win - 1)]);
  kl[0] = compute_kl(f1, f2);
  // i = 2
  f1 = density_cpp(x[Range(2 * jump, 2 * jump + win - 1)]);
  kl[1] = compute_kl(f2, f1);

  for (int i = 2; i < win_num - 1; i++) {

    if (Progress::check_abort())
      return -1.0;

    p.increment();

    if (kl_s < lambda) {
      f1 = density_cpp(x[Range(i * jump, i * jump + win - 1)]);
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2);
      kl_s = kl[i];
      lambda = lambda_p * (kl[i - 2] + eps);
    } else {
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2); // f1 = of normal
      kl_s = kl[i];
    }
  }

  return kl;
}




