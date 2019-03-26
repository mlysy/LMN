#' Profile loglikelihood for the LMN model.
#'
#' Calculate the loglikelihood of the LMN model defined in \code{\link{lmn.suff}} at the MLE \code{Beta = Beta.hat} and \code{Sigma = Sigma.hat}.
#'
#' @template param-suff
#' @template param-noSigma
#'
#' @return Scalar -- the calculated value of the profile loglikelihood.
#' @example examples/lmn-prof.R
#' @export
lmn.prof <- function(suff, noSigma = FALSE) {
  if(class(suff) != "lmn_suff") {
    stop("suff must be an object of class 'lmn_suff'.")
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
