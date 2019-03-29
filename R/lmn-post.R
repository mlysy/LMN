#' Parameters of the posterior conditional distribution of an LMN model.
#'
#' Calculates the parameters of the LMN model's Matrix-Normal Inverse-Wishart (MNIW) conjugate posterior distribution  (see \strong{Details}).
#'
#' @template param-suff
#' @param prior A list with elements \code{Lambda}, \code{Omega}, \code{Psi}, \code{nu} as returned by \code{\link{lmn.prior}}.
#'
#' @return A list with elements named as in \code{prior} specifying the parameters of the posterior MNIW distribution.  Elements \code{Omega = NA} and \code{nu = NA} specify that parameters \code{Beta = 0} and \code{Sigma = diag(q)}, respectively, are known and not to be estimated.
#'
#' @template details-mniw
#' @details The posterior MNIW distribution is required to be a proper distribution, but the prior is not.  For example, \code{prior = NULL} corresponds to the noninformative prior
#' \deqn{
#' \pi(B, \boldsymbol{\Sigma}) \sim |\boldsymbol{Sigma}|^{-(q+1)/2}.
#' }{
#' \pi(B, \Sigma) ~ |\Sigma|^{-(q+1)/2}.
#' }
#' @example examples/lmn-post.R
#' @export
lmn.post <- function(suff, prior) {
  if(class(suff) != "lmn_suff") {
    stop("suff must be an object of class 'lmn_suff'.")
  }
  n <- suff$n
  Betahat <- suff$Bhat
  T <- suff$T
  S <- suff$S
  noBeta <- is.null(T)
  p <- ifelse(noBeta, 0, nrow(Betahat))
  q <- nrow(S)
  ## # get prior
  ## if(is.null(prior) ||
  ##    !all(sort(names(prior)) == .PriorNames) ||
  ##    any(sapply(prior, is.null)) || prior$Omega[1] == 0) {
  ##   # try to avoid passing prior by reference to .DefaultPrior
  ##   prior <- .DefaultPrior(prior, p, q, noSigma)
  ## }
  Lambda <- prior$Lambda
  Omega <- prior$Omega
  Psi <- prior$Psi
  nu <- prior$nu
  # change noBeta to make consistent with prior
  if(!noBeta) noBeta <- (length(Omega) == 1) && is.na(Omega)
  noSigma <- is.na(nu)
  ## if(!noSigma) noSigma <- is.na(nu)
  if(noBeta && noSigma) {
    ## warning("Beta and Sigma both known.  calc.prior ignored.")
    return(list(Lambda = 0, Omega = NA, Psi = 0, nu = NA))
  }
  # calculate posterior parameters
  if(noBeta) {
    Lambdahat <- 0
    ## Lambdahat <- NA
    Omegahat <- NA
    noOmega <- TRUE
  } else {
    Omegahat <- T+Omega
    Lambdahat <- solveV(Omegahat, T%*%Betahat + Omega%*%Lambda)
    noOmega <- all(Omega == 0)
  }
  if(noSigma) {
    Psihat <- S
    ## Psihat <- 0
    nuhat <- NA
  } else {
    Psihat <- Psi + S
    nuhat <- nu + (n-p*noOmega)
  }
  if(!noBeta) {
    Psihat <- Psihat + crossprod(Betahat, T%*%Betahat) +
      crossprod(Lambda,Omega%*%Lambda) -
        crossprod(Lambdahat,Omegahat%*%Lambdahat)
  }
  list(Lambda = Lambdahat, Omega = Omegahat, Psi = Psihat, nu = nuhat)
  ## out <- list(Lambda = Lambdahat, Omega = Omegahat,
  ##             Psi = Psihat, nu = nuhat)
  ## if(calc.prior) {
  ##   prior <- list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
  ##   out <- c(out, list(prior = prior))
  ## }
  ## out
}

#--- scratch -------------------------------------------------------------------

## if(FALSE) {
##   Lambda <- prior$Lambda
##   Omega <- prior$Omega
##   nu <- prior$nu
##   Psi <- prior$Psi
##   if(is.null(nu)) nu <- 0
##   if(is.null(Psi)) Psi <- 0
##   if(!noBeta) {
##     noBeta <- (length(Omega) == 1) && is.na(Omega)
##   }
##   if(!noSigma) {
##     noSigma <- is.na(nu)
##   }
##   if(noBeta && noSigma) {
##     warning("Beta and Sigma both known.  calc.prior ignored.")
##     return(list(Lambda = 0, Omega = NA, Psi = 0, nu = NA))
##   }
##   if(noBeta) {
##     Lambda <- 0
##     Omega <- NA
##     noOmega <- TRUE
##   } else {
##     if(is.null(Lambda)) Lambda <- matrix(0,p,q)
##     Omega <- prior$Omega
##     if(is.null(Omega)) Omega <- 0
##     noOmega <- all(Omega == 0)
##     if(noOmega) Omega <- matrix(0,p,p)
##   }
##   if(noSigma) {
##     nu <- NA
##   }
## }
