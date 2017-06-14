#' Profile Log-Likelihood for the LMN Model
#'
#' @details The model is defined as
#' \deqn{Y ~ MNorm(X*Beta, V, Sigma),}
#' where MNorm is the Matrix-Normal distribution, i.e.
#' \deqn{vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).}
#' @param suff result of call to lmn.suff (i.e. avoid calculations here)
#' @param noSigma logical. if true assumes that \code{Sigma = diag(ncol(Y))}.
#' @return the calculated profile log-likelihood.
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
