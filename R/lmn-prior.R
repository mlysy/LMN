#' Conjugate prior specification for LMN models.
#'
#' The conjugate prior for LMN models is the Matrix-Normal Inverse-Wishart (MNIW) distribution.  This convenience function converts a partial MNIW prior specification into a full one.
#'
#' @param p Integer specifying row dimension of \code{Beta}.  \code{p = 0} corresponds to no \code{Beta} in the model, i.e., \code{X = 0} in \code{\link{lmn.suff}}.
#' @param q Integer specifying the dimension of \code{Sigma}.
#' @param Lambda Mean parameter for \code{Beta}.  Either:
#' \itemize{
#'   \item A \code{p x q} matrix.
#'   \item A scalar, in which case \code{Lambda = matrix(Lambda, p, q)}.
#'   \item Missing, in which case \code{Lambda = matrix(0, p, q)}.
#' }
#' @param Omega Row-wise precision parameter for \code{Beta}.  Either:
#' \itemize{
#'   \item A \code{p x p} matrix.
#'   \item A scalar, in which case \code{Omega = diag(rep(Omega,p))}.
#'   \item Missing, in which case \code{Omega = matrix(0, p, p)}.
#'   \item \code{NA}, which signifies that \code{Beta} is known, in which case the prior is purely Inverse-Wishart on \code{Sigma} (see \strong{Details}).
#' }
#' @param Psi Scale parameter for \code{Sigma}.  Either:
#' \itemize{
#'   \item A \code{q x q} matrix.
#'   \item A scalar, in which case \code{Psi = diag(rep(Psi,q))}.
#'   \item Missing, in which case \code{Psi = matrix(0, q, q)}.
#' }
#' @param nu Degrees-of-freedom parameter for \code{Sigma}.  Either a scalar, missing (defaults to \code{nu = 0}), or \code{NA}, which signifies that \code{Sigma = diag(q)} is known, in which case the prior is purely Matrix-Normal on \code{Beta} (see \strong{Details}).
#'
#' @template details-mniw
#'
#' @return A list with elements \code{Lambda}, \code{Omega}, \code{Psi}, \code{nu} with the proper dimensions specified above, except possibly \code{Omega = NA} or \code{nu = NA} (see \strong{Details}).
#'
#' @example examples/lmn-prior.R
#'
#' @export
lmn.prior <- function(p, q, Lambda, Omega, Psi, nu) {
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
