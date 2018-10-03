#' Log-likelihood function for LMN models.
#'
#' @param Beta an (p x q) matrix (see details).
#' @param Sigma an (q x q) matrix (see details).
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
#' 
#' @details The model is defined as
#' \deqn{Y ~ MNorm(X*Beta, V, Sigma),}
#' where MNorm is the Matrix-Normal distribution, i.e.
#' \deqn{vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).}
#' 
#' @return The calculated log-likelihood.
#' 
#' @examples
#' ## Data
#' Y = matrix(rnorm(100),50,2)
#' 
#' ## Scalar V example
#' X = matrix(1,50,1)
#' q = ncol(X)
#' V = 0.5
#' suff = lmn.suff(Y, X, V=V)
#' Beta = matrix(c(1,0.5),1,2)
#' Sigma = diag(2)
#' lmn.loglik(Beta, Sigma, suff, Y, X, V)
#' 
#' @export
lmn.loglik <- function(Beta, Sigma, suff, Y, X, V, acf) {
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, acf = acf)
  }
  # get sufficient statistics for likelihood evaluation
  n <- suff$n
  S <- suff$S
  q <- nrow(S)
  ldV <- suff$ldV
  Beta.hat <- suff$Beta.hat
  T <- suff$T
  noBeta <- is.null(Beta.hat)
  # log-likelihood calculation
  if(!noBeta) {
    Z <- Beta-Beta.hat
    S <- S + crossprod(Z, T %*% Z)
  }
  -0.5 * (sum(diag(solve(Sigma,S))) + q*ldV + n*ldet(Sigma) + n*q * log(2*pi))
}
