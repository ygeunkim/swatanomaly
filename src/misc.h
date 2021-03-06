#ifndef MISC_H
#define MISC_H

double sum_sq(NumericVector x);

NumericMatrix row_erase(NumericMatrix x, IntegerVector rowID);

IntegerVector seq_rcpp(int from, int to);

NumericMatrix sub_mat(NumericMatrix x, IntegerVector row, IntegerVector col);

LogicalVector rep_bool(bool x, int n);

NumericMatrix rbind_mat(NumericMatrix x, NumericMatrix y);

#endif
