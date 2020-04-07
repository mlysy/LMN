#' Conjugate prior specification for LMN models.
#'
#' The conjugate prior for LMN models is the Matrix-Normal Inverse-Wishart (MNIW) distribution.  This convenience function converts a partial MNIW prior specification into a full one.
#'
#' @param p Integer specifying row dimension of `Beta`.  `p = 0` corresponds to no `Beta` in the model, i.e., `X = 0` in [lmn_suff()].
#' @param q Integer specifying the dimension of `Sigma`.
#' @param Lambda Mean parameter for `Beta`.  Either:
#'
#' - A `p x q` matrix.
#' - A scalar, in which case `Lambda = matrix(Lambda, p, q)`.
#' - Missing, in which case `Lambda = matrix(0, p, q)`.
#'
#' @param Omega Row-wise precision parameter for `Beta`.  Either:
#'
#' - A `p x p` matrix.
#' - A scalar, in which case `Omega = diag(rep(Omega,p))`.
#' - Missing, in which case `Omega = matrix(0, p, p)`.
#' - `NA`, which signifies that `Beta` is known, in which case the prior is purely Inverse-Wishart on `Sigma` (see **Details**).
#'
#' @param Psi Scale parameter for `Sigma`.  Either:
#'
#' - A `q x q` matrix.
#' - A scalar, in which case `Psi = diag(rep(Psi,q))`.
#' - Missing, in which case `Psi = matrix(0, q, q)`.
#'
#' @param nu Degrees-of-freedom parameter for `Sigma`.  Either a scalar, missing (defaults to `nu = 0`), or `NA`, which signifies that `Sigma = diag(q)` is known, in which case the prior is purely Matrix-Normal on `Beta` (see **Details**).
#'
#' @template details-mniw
#'
#' @return A list with elements `Lambda`, `Omega`, `Psi`, `nu` with the proper dimensions specified above, except possibly `Omega = NA` or `nu = NA` (see **Details**).
#'
#' @example examples/lmn_prior.R
#'
#' @export
lmn_prior <- function(p, q, Lambda, Omega, Psi, nu) {
  noBeta <- p == 0
  # assign defaults
  if(missing(Lambda)) Lambda <- 0
  if(missing(Omega)) Omega <- 0
  if(missing(Psi)) Psi <- 0
  if(missing(nu)) nu <- 0
  if(!noBeta) {
    noBeta <- (length(Omega) == 1) && is.na(Omega)
  }
  if(noBeta) {
    Lambda <- 0
    Omega <- NA
  } else {
    if(length(Lambda) == 1) {
      Lambda <- matrix(Lambda[1],p,q)
    } else if(!.validDims(Lambda, c(p,q))) {
      stop("Lambda has invalid dimensions.")
    }
    if(length(Omega) == 1) {
      Omega <- diag(Omega[1],p)
    } else if(!.validDims(Omega, c(p,p))) {
      stop("Omega has invalid dimensions.")
    }
  }
  noSigma <- (length(nu) == 1) && is.na(nu)
  if(noSigma) {
    Psi <- 0
    nu <- NA
  } else {
    if(length(Psi) == 1) {
      Psi <- diag(Psi[1], q)
    } else if(!.validDims(Psi, c(q,q))){
      stop("Psi has invalid dimensions.")
    }
    if(length(nu) != 1) stop("nu has invalid dimensions.")
  }
  list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
}

.validDims <- function(X, dm) {
  isTRUE(all.equal(dim(X), dm, check.attributes = FALSE))
}
