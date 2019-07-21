#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;
#include "misc.h"

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

// [[Rcpp::export]]
int find_support(double x1, NumericVector x2) {
  int x_supp = 0;

  for (int i = 0; i < x2.size(); i++) {
    if (x2[i] > x1)
      break;
    x_supp++;
  }

  return x_supp;
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
  int x1_supp = 0;

  for (int i = 0; i < f1.nrow() - 1; i++) {

    x1_supp = find_support(f1(i, 0), f2(_, 1));

    if (x1_supp == 0) {
      sum_kl = 0;
    } else {
      px = f1(i, 1) * abs(f1(i + 1, 0) - f1(i, 0));
      qx = f2(x1_supp - 1, 1) * abs(f1(i + 1, 0) - f1(i, 0));

      if (px <= qx) {

        sum_kl = 0;

      } else {

        if (px == 0 | qx == 0) {
          sum_kl = 0;
        } else {
          sum_kl = px * log(px / qx);
        }
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
//' @return List,
//' First element is kl divergence named divergence.
//' Second element is threshold (lambda) for detecting anomaly named threshold.
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
List kl_dynamic(NumericVector x, int win, int jump, double lambda_p, double eps, bool display_progress = false) {
  int n = x.size();
  int win_num = (n - win) / jump + 1;
  NumericMatrix f1(512, 2);
  NumericMatrix f2(512, 2);

  // double lambda = lambda_p * eps;
  NumericVector lambda = rep(lambda_p * eps, win_num - 1);
  double kl_s = 0;

  Progress p(win_num - 3, display_progress);

  NumericVector kl(win_num - 1);
  LogicalVector anomaly(win_num - 1);

  // i = 1
  f1 = density_cpp(x[Range(0, win - 1)]);
  f2 = density_cpp(x[Range(jump, jump + win - 1)]);
  kl[0] = compute_kl(f1, f2);
  anomaly[0] = kl[0] > lambda[0];
  // i = 2
  f1 = density_cpp(x[Range(2 * jump, 2 * jump + win - 1)]);
  kl[1] = compute_kl(f2, f1);
  anomaly[1] = kl[1] > lambda[0];

  for (int i = 2; i < win_num - 1; i++) {

    if (Progress::check_abort())
      return -1.0;

    p.increment();

    if (kl_s < lambda[i]) {
      f1 = density_cpp(x[Range(i * jump, i * jump + win - 1)]);
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2);
      kl_s = kl[i];
      lambda[Range(i + 1, win_num - 1)] = rep(lambda_p * (kl[i - 2] + eps), win_num - 1 - i);

      anomaly[i] = kl[i] > lambda[i];
    } else {
      f2 = density_cpp(x[Range((i + 1) * jump, (i + 1) * jump + win - 1)]);

      kl[i] = compute_kl(f1, f2); // f1 = of normal
      kl_s = kl[i];

      anomaly[i] = kl[i] > lambda[i];
    }
  }

  List kl_alg = List::create(Named("divergence") = kl, Named("threshold") = lambda, Named("anomaly") = anomaly);

  return kl_alg;
}

//' Matching KL divergence label to individual observation
//'
//' @description
//' Give KL divergence anomaly prediction of each window to individual observation.
//' @param d LogicalVector anomaly of \code{\link{kl_dynamic}} or the result of detection by fixed algorithm
//' @param win int window size.
//' @return NumericVector of number identical to the original series except the last window.
//' @details
//' This function is not appropriate when jump option is used.
//' In other words, use only when the series has been partitioned.
//' @seealso
//'     \code{\link{kl_fix}}
//'     \code{\link{kl_dynamic}}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
LogicalVector match_kl(LogicalVector d, int win) {
  int dlen = d.size();
  int n = dlen * (win - 1);
  LogicalVector x(n);

  for (int i = 0; i < dlen; i++) {
    // window loop
    for (int j = 0; j < win; j++) {
      // in window
      x[i * win + j] = d[i];
    }
  }

  return x;
}

//' Online KL Algorithm
//'
//' @description
//' This function implement Dynamic lambda algorithm by online.
//' @param x NumericVector. univariate data set that consits of non-anomaly.
//' @param newx NumericVector. univariate data set that has possibility of anomaly.
//' @param win int window size.
//' @param jump int jump size for sliding window.
//' @param lambda_p double initializing lambda_p for the threshold.
//' @param eps double initializing epsilon for the threshold.
//' @param display_progress If TRUE, display a progress bar. By default, FALSE.
//' @return List,
//' First element is kl divergence named divergence.
//' Second element is threshold (lambda) for detecting anomaly named threshold.
//' @details
//' This is an online version for dynamic lambda algorithm.
//' In this setting, normal data set is given. We keep updating new data set that has possibility of anomaly.
//' This function tries to detect anomaly in this updated set.
//' First, estimate kernel in the normal set.
//' Next, estimate kernel in the first window of new data set.
//' Compute the KL and see if the window is anomaly.
//' If it is normal, estimate kernel in the normal set including the first window.
//' Otherwise, just keep the former kernel.
//' Compute the KL from second window and see if this window is anomaly.
//' Repeat this procedure while updating lambda.
//' @seealso
//'    \link[stats]{density.default}
//'     \code{\link{est_density}}
//'     \code{\link{density_cpp}}
//'     \code{\link{compute_kl}}
//'     \code{\link{kl_dynamic}}
//' @references Cho, J., Tariq, S., Lee, S., Kim, Y. G., & Woo, S. (2019). \emph{Contextual Anomaly Detection by Correlated Probability Distributions using Kullback-Leibler Divergence}. Workshop on Mining and Learning From Time Series. \url{http://doi.org/10.1145/nnnnnnn.nnnnnnn}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List kl_online(
  NumericVector x,
  NumericVector newx,
  int win,
  int jump,
  double lambda_p,
  double eps,
  bool display_progress
) {
  int n = x.size();
  int n_new = newx.size();
  int win_num = (n - win) / jump + 1;
  int win_new = (n_new - win) / jump + 1;

  NumericMatrix f1(512, 2);
  NumericMatrix f2(512, 2);

  // double lambda = lambda_p * eps;
  NumericVector lambda = rep(lambda_p * eps, win_new);
  // double kl_s = 0;

  Progress p(win_num - 1, display_progress);

  NumericVector kl(win_new);
  LogicalVector anomaly(win_new);
  NumericVector normal_win(n + n_new);
  normal_win[Range(0, n - 1)] = x;
  normal_win[Range(n, n + n_new - 1)] = newx;
  int test = 0;

  // i = 1
  f1 = density_cpp(normal_win[Range(0, n - 1)]);
  f2 = density_cpp(newx[Range(0, win - 1)]);
  kl[0] = compute_kl(f1, f2);
  anomaly[0] = kl[0] > lambda[0];

  for(int i = 1; i < win_new; i++) {

    if (Progress::check_abort())
      return -1.0;

    p.increment();

    if (!anomaly[i - 1]) {
      normal_win[Range(n + test * win, n + test * win + win - 1)] = newx[Range((i - 1) * jump, (i - 1) * jump + win - 1)];
      f1 = density_cpp(normal_win[Range(0, n + test * win + win - 1)]);
      test++;
    }

    if (kl[i - 1] < lambda[i - 1])
      lambda[Range(i, win_new)] = rep(lambda_p * (kl[i - 2] + eps), win_new - i + 1);

    f2 = density_cpp(newx[Range(i * jump, i * jump + win - 1)]);
    kl[i] = compute_kl(f1, f2);
    anomaly[i] = kl[i] > lambda[i];
  }

  List kl_alg = List::create(Named("divergence") = kl, Named("threshold") = lambda, Named("anomaly") = anomaly);

  return kl_alg;
}


