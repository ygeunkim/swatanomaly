#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;
#include "misc.h"

//' Sliding window for NND
//'
//' @description
//' For the chosen partition, get NND for every pair of window.
//' This function is constructed for the other function.
//' @param data NumericMatrix. data to be calculated NND.
//' @param win int. window size.
//' @param jump int. shift size.
//' @param partition NumericVector. indices that can indicate partitions.
//' @param base_id int. An index of chosen partition among partiton vector.
//' @param d Function. distance function.
//' @return NumericVector. NND vector for each window in the chosen partition.
//' @seealso
//'     \code{\link{nnd_normal}}
//'     \code{\link{pred_nnd}}
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
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
  Function d,
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
//' @details
//' In the window in new data set, scan the original set sliding window.
//' Minimum distance is the NND of that window.
//' Repeat for every window.
//' @useDynLib swatanomaly
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
NumericVector pred_nnd(
  NumericMatrix data,
  NumericMatrix newdata,
  int win,
  int jump,
  Function d,
  bool display_progress = false
) {
  int n = data.nrow();
  int new_win = (newdata.nrow() - win) / jump + 1;
  NumericMatrix dat_new(n + win, data.ncol());
  NumericVector nnd_new(new_win);

  IntegerVector id(n + win);
  id[Range(0, n - 1)] = rep(0, n);
  id[Range(n, n + win - 1)] = rep(1, win);

  for (int i = 0; i < new_win; i++) {
    dat_new = rbind_mat(data, newdata(Range(i * jump, i * jump + win - 1), _));
    nnd_new[i] = as<double>(partnnd(dat_new, win, jump, id, rep(1, win), d));
  }

  return nnd_new;
}

