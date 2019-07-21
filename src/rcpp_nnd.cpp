#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;
#include <iostream>
#include "misc.h"
#include "distnnd.h"

//' Sliding window for NND
//'
//' @description
//' For the chosen partition, get NND for every pair of window.
//' This function is constructed for the other function.
//' @param data NumericMatrix. data to be calculated NND.
//' @param win int. window size.
//' @param jump int. shift size.
//' @param from int. First index of chosen partition.
//' @param to int. Last index of chosen partition.
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
  int from,
  int to
) {
  IntegerVector id = seq_len(data.nrow()) - 1;
  IntegerVector x_id = id[id >= from & id <= to];
  IntegerVector y_id = id[id < from | id > to];

  int px = data.ncol();
  NumericMatrix x(x_id.size(), px);
  NumericMatrix y(y_id.size(), px);
  x = sub_mat(data, x_id, seq_len(px) - 1);
  y = sub_mat(data, y_id, seq_len(px) - 1);

  int x_win = (x_id.size() - win) / jump + 1;
  NumericVector x_nnd(x_win);
  int y_win = (y_id.size() - win) / jump + 1;
  NumericVector candid(y_win);

  for (int i = 0; i < x_win; i++) {
    for (int j = 0; j < y_win; j++) {
      candid[j] = compute_euc(sub_mat(x, seq(i * jump, i * jump + win - 1), seq_len(px) - 1),
                             sub_mat(y, seq(j * jump, j * jump + win - 1), seq_len(px) - 1));
    }
    x_nnd[i] = min(candid);
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
  bool display_progress = false
) {
  int n = data.nrow();
  int part_num = n / part; // obs number in each partition
  int part_rem = n % part; // remainder in last partition

  Progress p(part, display_progress);

  int start = 0;
  int end = part_num - 1;

  int win_num = (part_num - win) / jump + 1; // window in partition
  int win_last = (part_num + part_rem - win) / jump + 1; // window in last partition

  NumericVector nnd(win_num * (part - 1) + win_last);

  for (int i = 0; i < (part - 1); i++) {
    if (Progress::check_abort())
      return -1.0;

    p.increment();
    nnd[Range(start, end)] = partnnd(data, win, jump, start, end);
    start = end + 1;
    end += part_num;
  }

  p.increment();
  end += part_rem;
  nnd[Range(start, end)] = partnnd(data, win, jump, start, end);

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
  bool display_progress = false
) {
  int n = data.nrow();
  int new_win = (newdata.nrow() - win) / jump + 1; // window number in newdata
  NumericMatrix dat_new(n + win, data.ncol()); // merge data + one window in newdata
  NumericVector nnd_new(new_win);

  Progress p(new_win, display_progress);

  int start = n;
  int end = n + win - 1;

  for (int i = 0; i < new_win; i++) {
    if (Progress::check_abort())
      return -1.0;

    p.increment();

    dat_new = rbind_mat(data,
                        sub_mat(newdata, seq(i * jump, i * jump + win - 1), seq_len(newdata.ncol()) - 1));
    nnd_new[i] = as<double>(partnnd(dat_new, win, jump, start, end));
  }

  return nnd_new;
}

