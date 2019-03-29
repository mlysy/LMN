#' Loglikelihood function for LMN models.
#'
#' @param Beta A \code{p x q} matrix of regression coefficients (see \code{\link{lmn.suff}}).
#' @param Sigma A \code{q x q} matrix of error variances (see \code{\link{lmn.suff}}).
#' @template param-suff
#' @return Scalar -- the value of the loglikelihood.
#' @example examples/lmn-loglik.R
#' @export
lmn.loglik <- function(Beta, Sigma, suff) {
  if(class(suff) != "lmn_suff") {
    stop("suff must be an object of class 'lmn_suff'.")
  }
  # get sufficient statistics for likelihood evaluation
  n <- suff$n
  S <- suff$S
  q <- nrow(S)
  ldV <- suff$ldV
  Bhat <- suff$Bhat
  T <- suff$T
  noBeta <- is.null(Bhat)
  # log-likelihood calculation
  if(!noBeta) {
    Z <- Beta-Bhat
    S <- S + crossprod(Z, T %*% Z)
  }
  -0.5 * (sum(diag(solveV(Sigma,S))) + q*ldV + n*ldet(Sigma) + n*q * log(2*pi))
}
