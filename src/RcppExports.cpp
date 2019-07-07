// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

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
// euc_dist
double euc_dist(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _swatanomaly_euc_dist(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(euc_dist(x, y));
    return rcpp_result_gen;
END_RCPP
}
// euc_nnd
NumericVector euc_nnd(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _swatanomaly_euc_nnd(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(euc_nnd(x, y));
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
// euc_pdf
NumericVector euc_pdf(NumericMatrix x, int partition, bool display_progress);
RcppExport SEXP _swatanomaly_euc_pdf(SEXP xSEXP, SEXP partitionSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type partition(partitionSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(euc_pdf(x, partition, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// nns_cpp
NumericVector nns_cpp(NumericMatrix data, int win, bool display_progress);
RcppExport SEXP _swatanomaly_nns_cpp(SEXP dataSEXP, SEXP winSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(nns_cpp(data, win, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// detect_nnd
LogicalVector detect_nnd(NumericMatrix data, int win, double thr);
RcppExport SEXP _swatanomaly_detect_nnd(SEXP dataSEXP, SEXP winSEXP, SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type win(winSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(detect_nnd(data, win, thr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_swatanomaly_sum_sq", (DL_FUNC) &_swatanomaly_sum_sq, 1},
    {"_swatanomaly_euc_dist", (DL_FUNC) &_swatanomaly_euc_dist, 2},
    {"_swatanomaly_euc_nnd", (DL_FUNC) &_swatanomaly_euc_nnd, 2},
    {"_swatanomaly_row_erase", (DL_FUNC) &_swatanomaly_row_erase, 2},
    {"_swatanomaly_seq_rcpp", (DL_FUNC) &_swatanomaly_seq_rcpp, 2},
    {"_swatanomaly_euc_pdf", (DL_FUNC) &_swatanomaly_euc_pdf, 3},
    {"_swatanomaly_nns_cpp", (DL_FUNC) &_swatanomaly_nns_cpp, 3},
    {"_swatanomaly_detect_nnd", (DL_FUNC) &_swatanomaly_detect_nnd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_swatanomaly(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
