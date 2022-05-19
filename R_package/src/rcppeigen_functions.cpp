// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// Output a matrix equal to the matrix product X*transpose(X) for matrix X
// [[Rcpp::export]]
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X){
  const int n = X.cols();
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}

// Output a matrix equal to the matrix product X*Y for matrices X and Y
// [[Rcpp::export]]
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
  return Eigen::MatrixXd(X*Y);
}

// Output a vector equal to the matrix product X*y for matrix X and vector y
// [[Rcpp::export]]
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y){
  return Eigen::VectorXd(X*y);
}
