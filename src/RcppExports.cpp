// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_euc
double compute_euc(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _swatanomaly_compute_euc(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_euc(x, y));
    return rcpp_result_gen;
END_RCPP
}
// detect
LogicalVector detect(NumericVector y, int win, int jump, double thr);
RcppExport SEXP _swatanomaly_detect(SEXP ySEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(detect(y, win, jump, thr));
    return rcpp_result_gen;
END_RCPP
}
// expand_label
LogicalVector expand_label(LogicalVector x, int win);
RcppExport SEXP _swatanomaly_expand_label(SEXP xSEXP, SEXP winSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    rcpp_result_gen = Rcpp::wrap(expand_label(x, win));
    return rcpp_result_gen;
END_RCPP
}
// aggregate_mts
NumericVector aggregate_mts(NumericMatrix x);
RcppExport SEXP _swatanomaly_aggregate_mts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregate_mts(x));
    return rcpp_result_gen;
END_RCPP
}
// density_cpp
NumericMatrix density_cpp(NumericVector x);
RcppExport SEXP _swatanomaly_density_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(density_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// find_support
int find_support(double x1, NumericVector x2);
RcppExport SEXP _swatanomaly_find_support(SEXP x1SEXP, SEXP x2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    rcpp_result_gen = Rcpp::wrap(find_support(x1, x2));
    return rcpp_result_gen;
END_RCPP
}
// compute_kl
double compute_kl(NumericMatrix f1, NumericMatrix f2);
RcppExport SEXP _swatanomaly_compute_kl(SEXP f1SEXP, SEXP f2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type f2(f2SEXP);
    rcpp_result_gen = Rcpp::wrap(compute_kl(f1, f2));
    return rcpp_result_gen;
END_RCPP
}
// kl_fix
NumericVector kl_fix(NumericVector x, int win, int jump, double lambda, bool display_progress);
RcppExport SEXP _swatanomaly_kl_fix(SEXP xSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP lambdaSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(kl_fix(x, win, jump, lambda, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// kl_dynamic
List kl_dynamic(NumericVector x, int win, int jump, double lambda_p, double eps, bool display_progress);
RcppExport SEXP _swatanomaly_kl_dynamic(SEXP xSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP lambda_pSEXP, SEXP epsSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_p(lambda_pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(kl_dynamic(x, win, jump, lambda_p, eps, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// match_kl
LogicalVector match_kl(LogicalVector d, int win);
RcppExport SEXP _swatanomaly_match_kl(SEXP dSEXP, SEXP winSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    rcpp_result_gen = Rcpp::wrap(match_kl(d, win));
    return rcpp_result_gen;
END_RCPP
}
// kl_online
List kl_online(NumericVector x, NumericVector newx, int win, int jump, double lambda_p, double eps, bool display_progress);
RcppExport SEXP _swatanomaly_kl_online(SEXP xSEXP, SEXP newxSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP lambda_pSEXP, SEXP epsSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type newx(newxSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_p(lambda_pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(kl_online(x, newx, win, jump, lambda_p, eps, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// sum_sq
double sum_sq(NumericVector x);
RcppExport SEXP _swatanomaly_sum_sq(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_sq(x));
    return rcpp_result_gen;
END_RCPP
}
// row_erase
NumericMatrix row_erase(NumericMatrix x, IntegerVector rowID);
RcppExport SEXP _swatanomaly_row_erase(SEXP xSEXP, SEXP rowIDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowID(rowIDSEXP);
    rcpp_result_gen = Rcpp::wrap(row_erase(x, rowID));
    return rcpp_result_gen;
END_RCPP
}
// seq_rcpp
IntegerVector seq_rcpp(int from, int to);
RcppExport SEXP _swatanomaly_seq_rcpp(SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_rcpp(from, to));
    return rcpp_result_gen;
END_RCPP
}
// sub_mat
NumericMatrix sub_mat(NumericMatrix x, IntegerVector row, IntegerVector col);
RcppExport SEXP _swatanomaly_sub_mat(SEXP xSEXP, SEXP rowSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(sub_mat(x, row, col));
    return rcpp_result_gen;
END_RCPP
}
// rep_bool
LogicalVector rep_bool(bool x, int n);
RcppExport SEXP _swatanomaly_rep_bool(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_bool(x, n));
    return rcpp_result_gen;
END_RCPP
}
// rbind_mat
NumericMatrix rbind_mat(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _swatanomaly_rbind_mat(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(rbind_mat(x, y));
    return rcpp_result_gen;
END_RCPP
}
// concat_vec
NumericVector concat_vec(NumericVector x, NumericVector y);
RcppExport SEXP _swatanomaly_concat_vec(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(concat_vec(x, y));
    return rcpp_result_gen;
END_RCPP
}
// partnnd
NumericVector partnnd(NumericMatrix data, int win, int jump, int from, int to);
RcppExport SEXP _swatanomaly_partnnd(SEXP dataSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(partnnd(data, win, jump, from, to));
    return rcpp_result_gen;
END_RCPP
}
// nnd_normal
NumericVector nnd_normal(NumericMatrix data, int part, int win, int jump, bool display_progress);
RcppExport SEXP _swatanomaly_nnd_normal(SEXP dataSEXP, SEXP partSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type part(partSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(nnd_normal(data, part, win, jump, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// pred_nnd
NumericVector pred_nnd(NumericMatrix data, NumericMatrix newdata, int win, int jump, bool display_progress);
RcppExport SEXP _swatanomaly_pred_nnd(SEXP dataSEXP, SEXP newdataSEXP, SEXP winSEXP, SEXP jumpSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type newdata(newdataSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< int >::type jump(jumpSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_nnd(data, newdata, win, jump, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_swatanomaly_compute_euc", (DL_FUNC) &_swatanomaly_compute_euc, 2},
    {"_swatanomaly_detect", (DL_FUNC) &_swatanomaly_detect, 4},
    {"_swatanomaly_expand_label", (DL_FUNC) &_swatanomaly_expand_label, 2},
    {"_swatanomaly_aggregate_mts", (DL_FUNC) &_swatanomaly_aggregate_mts, 1},
    {"_swatanomaly_density_cpp", (DL_FUNC) &_swatanomaly_density_cpp, 1},
    {"_swatanomaly_find_support", (DL_FUNC) &_swatanomaly_find_support, 2},
    {"_swatanomaly_compute_kl", (DL_FUNC) &_swatanomaly_compute_kl, 2},
    {"_swatanomaly_kl_fix", (DL_FUNC) &_swatanomaly_kl_fix, 5},
    {"_swatanomaly_kl_dynamic", (DL_FUNC) &_swatanomaly_kl_dynamic, 6},
    {"_swatanomaly_match_kl", (DL_FUNC) &_swatanomaly_match_kl, 2},
    {"_swatanomaly_kl_online", (DL_FUNC) &_swatanomaly_kl_online, 7},
    {"_swatanomaly_sum_sq", (DL_FUNC) &_swatanomaly_sum_sq, 1},
    {"_swatanomaly_row_erase", (DL_FUNC) &_swatanomaly_row_erase, 2},
    {"_swatanomaly_seq_rcpp", (DL_FUNC) &_swatanomaly_seq_rcpp, 2},
    {"_swatanomaly_sub_mat", (DL_FUNC) &_swatanomaly_sub_mat, 3},
    {"_swatanomaly_rep_bool", (DL_FUNC) &_swatanomaly_rep_bool, 2},
    {"_swatanomaly_rbind_mat", (DL_FUNC) &_swatanomaly_rbind_mat, 2},
    {"_swatanomaly_concat_vec", (DL_FUNC) &_swatanomaly_concat_vec, 2},
    {"_swatanomaly_partnnd", (DL_FUNC) &_swatanomaly_partnnd, 5},
    {"_swatanomaly_nnd_normal", (DL_FUNC) &_swatanomaly_nnd_normal, 5},
    {"_swatanomaly_pred_nnd", (DL_FUNC) &_swatanomaly_pred_nnd, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_swatanomaly(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
