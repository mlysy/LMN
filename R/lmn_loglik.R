#' Loglikelihood function for LMN models.
#'
#' @param Beta A `p x q` matrix of regression coefficients (see [lmn_suff()]).
#' @param Sigma A `q x q` matrix of error variances (see [lmn_suff()]).
#' @template param-suff
#' @return Scalar; the value of the loglikelihood.
#' @example examples/lmn_loglik.R
#' @export
lmn_loglik <- function(Beta, Sigma, suff) {
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
  IP <- solveV(Sigma, S, ldV = TRUE)
  -0.5 * (sum(diag(IP$y)) + q*ldV + n*IP$ldV + n*q * log(2*pi))
  ## -0.5 * (sum(diag(solveV(Sigma,S))) + q*ldV + n*ldet(Sigma) + n*q * log(2*pi))
}
