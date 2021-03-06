// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _ScaleSpikeSlab_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// fcprd
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X);
RcppExport SEXP _ScaleSpikeSlab_fcprd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fcprd(X));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prod
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y);
RcppExport SEXP _ScaleSpikeSlab_cpp_prod(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prod(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mat_vec_prod
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y);
RcppExport SEXP _ScaleSpikeSlab_cpp_mat_vec_prod(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mat_vec_prod(X, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ScaleSpikeSlab_rcpp_hello_world", (DL_FUNC) &_ScaleSpikeSlab_rcpp_hello_world, 0},
    {"_ScaleSpikeSlab_fcprd", (DL_FUNC) &_ScaleSpikeSlab_fcprd, 1},
    {"_ScaleSpikeSlab_cpp_prod", (DL_FUNC) &_ScaleSpikeSlab_cpp_prod, 2},
    {"_ScaleSpikeSlab_cpp_mat_vec_prod", (DL_FUNC) &_ScaleSpikeSlab_cpp_mat_vec_prod, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ScaleSpikeSlab(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
