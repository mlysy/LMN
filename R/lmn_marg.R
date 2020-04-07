#' Marginal log-posterior for the LMN model.
#'
#' @template param-suff
#' @param prior A list with elements `Lambda`, `Omega`, `Psi`, `nu` corresponding to the parameters of the prior MNIW distribution.  See [lmn_prior()].
#' @param post A list with elements `Lambda`, `Omega`, `Psi`, `nu` corresponding to the parameters of the posterior MNIW distribution.  See [lmn_post()].
#' @return The scalar value of the marginal log-posterior.
#'
#' @example examples/lmn_marg.R
#' @export
lmn_marg <- function(suff, prior, post) {
  if(class(suff) != "lmn_suff") {
    stop("suff must be an object of class 'lmn_suff'.")
  }
  n <- suff$n
  Betahat <- suff$Bhat
  S <- suff$S
  ldV <- suff$ldV
  noBeta <- is.null(suff$T)
  p <- ifelse(noBeta, 0, nrow(Betahat))
  q <- nrow(S)
  # posterior MNIW parameters
  Omegahat <- post$Omega
  Psihat <- post$Psi
  nuhat <- post$nu
  noBeta <- (length(Omegahat) == 1) && is.na(Omegahat)
  noSigma <- is.na(nuhat)
  # prior MNIW parameters
  nu <- prior$nu
  Psi <- prior$Psi
  Omega <- prior$Omega
  if(noBeta) {
    noOmega <- TRUE
  } else {
    noOmega <- all(Omega == 0)
  }
  if(noSigma) {
    noPsi <- TRUE
  } else {
    noPsi <- all(Psi == 0)
  }
  # marginal log-posterior
  lp <- 0
  if(!noBeta) {
    lp <- lp + ldet(Omegahat)
  }
  lp <- q * (lp + ldV)
  if(!noSigma) {
    lp <- -.5 * (nuhat * (ldet(Psihat) - q*log(2)) + lp)
    lp <- lp + lmgamma(.5*nuhat, q)
  } else {
    lp <- -.5 * (lp + sum(diag(Psihat)))
  }
  # normalize by prior
  lpi <- 0
  if(!noBeta && !noOmega) {
    lpi <- lpi + ldet(Omega)
  }
  lpi <- q * (lpi - (n - p*noOmega) * log(2*pi))
  if(!noSigma && !noPsi) {
    lpi <- -.5 * (nu * (ldet(Psi) - q*log(2)) + lpi)
    lpi <- lpi + lmgamma(.5*nu, q)
  } else {
    lpi <- -.5 * lpi
  }
  lp - lpi
}
