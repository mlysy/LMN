// Durbin Levinson algorithm for inner products of the form
// M = X' toeplitz(acf)^{-1} Y
//
// Two versions: with and without eigen

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "DurbinLevinson.h"

//[[Rcpp::export("DurbinLevinsonEigen")]]
Rcpp::List DurbinLevinson_Eigen(Eigen::MatrixXd X, Eigen::MatrixXd Y,
				     Eigen::VectorXd acf,
				     int calcMode = 1) {
  int n, d, k;
  n = acf.size();
  d = X.rows();
  k = calcMode == 1 ? X.rows() : Y.rows();
  // output
  MatrixXd M(calcMode != 2 ? d : 1, k);
  double ldV = 0.0;
  // tmp
  VectorXd phi(n);
  VectorXd phi2(n);
  VectorXd rx(d);
  VectorXd ry(k);

  DurbinLevisonEigen(M, ldV, X, Y, acf, phi, phi2, rx, ry, calcMode);
  return List::create(_["IP"] = wrap(M), _["ldV"] = wrap(ldV));
}

//[[Rcpp::export("DurbinLevinsonBase")]]
Rcpp::List DurbinLevinson_Base(NumericMatrix X, NumericMatrix Y,
			       NumericVector acf, int calcMode = 1) {
  int n, d, k;
  n = acf.length();
  d = X.ncol();
  k = calcMode == 1 ? d : Y.ncol();
  // output
  NumericMatrix M(calcMode != 2 ? d : 1, k);
  double ldV = 0.0;
  // tmp
  double *phi = new double[n];
  double *phi2 = new double[n];
  double *rx = new double[d];
  double *ry = new double[k];

  DurbinLevinsonBase(REAL(M), ldV, REAL(X), REAL(Y), REAL(acf),
		     phi, phi2, rx, ry, n, d, k, calcMode);
  delete [] phi;
  delete [] phi2;
  delete [] rx;
  delete [] ry;
  return List::create(_["IP"] = wrap(M), _["ldV"] = wrap(ldV));
}
