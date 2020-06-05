#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

/// Inner product with the inverse of a positive-definite symmetric matrix.
///
/// @param[in] V Symmetric positive-definite matrix of size `N x N`.
/// @param[in] Z Matrix of size `N x q`.
///
/// @return An `Rcpp::List` with named elements `IP = Z' V^{-1} Z` and `ldV = log(det(V))`.
// [[Rcpp::export]]
Rcpp::List CholeskyIP(Eigen::MatrixXd V, Eigen::MatrixXd Z) {
  int N = Z.rows();
  int q = Z.cols();
  LLT<MatrixXd> cholV(N);
  MatrixXd iCZ(N,q);
  MatrixXd IP(q,q);
  double ldV = 0.0;
  
  cholV.compute(V);
  iCZ = Z;
  cholV.matrixL().solveInPlace(iCZ);
  IP.noalias() = iCZ.transpose() * iCZ;
  for(int ii=0; ii<N; ii++) {
    ldV += log(cholV.matrixL()(ii,ii));
  }
  ldV *= 2.0;
  return List::create(_["IP"] = IP, _["ldV"] = ldV);
}
