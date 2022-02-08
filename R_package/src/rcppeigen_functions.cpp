// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// Fast crossproduct of single matrix
// [[Rcpp::export]]
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X){
    const int n = X.cols();
    return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}
// Fast crossproduct of two matrices
// [[Rcpp::export]]
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
  return Eigen::MatrixXd(X*Y);
}
// Fast product of a matrix and vector
// [[Rcpp::export]]
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y){
  return Eigen::VectorXd(X*y);
}

