#' Profile log-likelihood for the LMN model.
#'
#' @details The model is defined as
#' \deqn{Y ~ MNorm(X*Beta, V, Sigma),}
#' where MNorm is the Matrix-Normal distribution, i.e.
#' \deqn{vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).}
#' @param suff result of call to lmn.suff (optional, to avoid calculations).
#' @param Y an (n x q) matrix.
#' @param X an (n x p) matrix. May also be passed as
#' \itemize{
#'    \item \code{X = 0}: in which case there is no intercept,
#'    \item \code{X != 0}: in which case a scaled intercept X = X * matrix(1, n, 1)
#'    is assumed.
#' }
#' @param V either: (1) an \code{n x n} full matrix, (2) an vector of length \code{n} such that \code{V = diag(V)}, (3) a scalar, such that \code{V = V * diag(n)}.
#' @param acf a vector of length n such that V = toeplitz(acf).
#' @param noSigma logical. If true assumes that \code{Sigma = diag(ncol(Y))}.
#' @return The calculated profile log-likelihood.
#' 
#' @examples
#' ## Data
#' Y = matrix(rnorm(100),50,2)
#' 
#' ## example where lmn.suff if first called
#' X = matrix(1,50,1)
#' V = exp(-seq(1:nrow(Y)))
#' suff = lmn.suff(Y, X, V=V)
#' lmn.prof(suff, Y, X, V)
#' 
#' ## lmn.prof call without suff gives same result (calculates suff internally)
#' lmn.prof(Y=Y, X=X, V=V)
#' 
#' @export
lmn.prof <- function(suff, Y, X, V, acf, noSigma = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, acf = acf)
  }
  n <- suff$n
  S <- suff$S
  ldV <- suff$ldV
  q <- nrow(S)
  if(!noSigma) {
    ll <- n*q * (1 + log(2*pi)) + n*ldet(S/n) + q*ldV
  } else {
    ll <- n*q * log(2*pi) + sum(diag(S)) + q*ldV
  }
  -.5 * ll
}
