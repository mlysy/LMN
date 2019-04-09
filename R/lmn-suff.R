#' Calculate the sufficient statistics of an LMN model.
#'
#' @param Y An \code{n x q} matrix of responses.
#' @param X An \code{N x p} matrix of covariates, where \code{N = n + npred} (see \strong{Details}). May also be passed as:
#' \itemize{
#'   \item A scalar, in which case the one-column covariate matrix is \code{X = X * matrix(1, N, 1)}.
#'   \item \code{X = 0}, in which case the mean of \code{Y} is known to be zero, i.e., no regression coefficients are estimated.
#' }
#' @param V,Vtype The between-observation variance specification.  Currently the following options are supported:
#' \itemize{
#'   \item \code{Vtype = "full"}: \code{V} is an \code{N x N} symmetric positive-definite matrix.
#'   \item \code{Vtype = "diag"}: \code{V} is a vector of length \code{N} such that \code{V = diag(V)}.
#'   \item \code{Vtype = "scalar"}: \code{V} is a scalar such that \code{V = V * diag(N)}.
#'   \item \code{Vtype = "acf"}: \code{V} is either a vector of length \code{N}, or an object of class \code{\link[SuperGauss]{Toeplitz}} such that \code{V = toeplitz(V)}.
#' }
#' For \code{V} specified as a matrix or scalar, \code{Vtype} is deduced automatically and need not be specified.
#' @param npred A nonnegative integer.  If positive, calculates sufficient statistics to make predictions for new responses. See \strong{Details}.
#' @return An S3 object of type \code{lmn_suff}, consisting of a list with elements:
#' \describe{
#'   \item{\code{Bhat}}{The \eqn{p \times q}{p x q} matrix \eqn{\hat{\boldsymbol{B}} = (\boldsymbol{X}'\boldsymbol{V}^{-1}\boldsymbol{X})^{-1}\boldsymbol{X}'\boldsymbol{V}^{-1}\boldsymbol{Y}}{B_hat = (X'V^{-1}X)^{-1}X'V^{-1}Y}.}
#'   \item{\code{T}}{The \eqn{p \times p}{p x p} matrix \eqn{\boldsymbol{T} = \boldsymbol{X}'\boldsymbol{V}^{-1}\boldsymbol{X}}{T = X'V^{-1}X}.}
#'   \item{\code{S}}{The \eqn{q \times q}{q x q} matrix \eqn{\boldsymbol{S} = (\boldsymbol{Y} - \boldsymbol{X} \hat{\boldsymbol{B}})'\boldsymbol{V}^{-1}(\boldsymbol{Y} - \boldsymbol{X} \hat{\boldsymbol{B}})}{S = (Y-X B_hat)'V^{-1}(Y-X B_hat)}.}
#'   \item{\code{ldV}}{The scalar log-determinant of \code{V}.}
#'   \item{\code{n}, \code{p}, \code{q}}{The problem dimensions, namely \code{n = nrow(Y)}, \code{p = nrow(Beta)} (or \code{p = 0} if \code{X = 0}), and \code{q = ncol(Y)}.}
#' }
#' In addition, when \code{npred > 0} and with \eqn{\boldsymbol{x}}{x}, \eqn{\boldsymbol{w}}{w}, and \eqn{v} defined in \strong{Details}:
#' \describe{
#'   \item{\code{Ap}}{The \code{npred x q} matrix \eqn{\boldsymbol{A}_p = \boldsymbol{w}'\boldsymbol{V}^{-1}\boldsymbol{Y}}{A_p = w'V^{-1}Y}.}
#'   \item{\code{Xp}}{The \code{npred x p} matrix \eqn{\boldsymbol{X}_p = \boldsymbol{x} - \boldsymbol{w}\boldsymbol{V}^{-1}\boldsymbol{X}}{X_p = x - w'V^{-1}X}.}
#'   \item{\code{Vp}}{The scalar \eqn{V_p = v - \boldsymbol{w}\boldsymbol{V}^{-1}\boldsymbol{w}}{V_p = v - w'V^{-1}w}.}
#' }
#'
#' @details
#' The multi-response normal linear regression model is defined as
#' \deqn{
#' \boldsymbol{Y} \sim \textrm{Matrix-Normal}(\boldsymbol{X}\boldsymbol{B}, \boldsymbol{V}, \boldsymbol{\Sigma}),
#' }{
#' Y ~ Matrix-Normal(X B, V, \Sigma),
#' }
#' where \eqn{\boldsymbol{Y}_{n \times q}}{Y_(n x q)} is the response matrix, \eqn{\boldsymbol{X}_{n \times p}}{X_(n x p)} is the covariate matrix, \eqn{\boldsymbol{B}_{p \times q}}{B_(p x q)} is the coefficient matrix, \eqn{\boldsymbol{V}_{n \times n}}{V_(n x n)} and \eqn{\boldsymbol{\Sigma}_{q \times q}}{\Sigma_(q x q)} are the between-row and between-column variance matrices, and the Matrix-Normal distribution is defined by the multivariate normal distribution
#' \eqn{
#' \textrm{vec}(\boldsymbol{Y}) \sim \mathcal{N}(\textrm{vec}(\boldsymbol{X}\boldsymbol{B}), \boldsymbol{\Sigma} \otimes \boldsymbol{V}),
#' }{
#' vec(Y) ~ N( vec(X B), \Sigma \%x\% V ),
#' }
#' where \eqn{\textrm{vec}(\boldsymbol{Y})}{vec(Y)} is a vector of length \eqn{nq} stacking the columns of of \eqn{\boldsymbol{Y}}{Y}, and \eqn{\boldsymbol{\Sigma} \otimes \boldsymbol{V}}{\Sigma \%x\% V} is the Kronecker product.
#'
#' The function \code{lmn.suff} returns everything needed to efficiently calculate the likelihood function
#' \deqn{\mathcal{L}(\boldsymbol{B}, \boldsymbol{\Sigma} \mid \boldsymbol{Y}, \boldsymbol{X}, \boldsymbol{V}) = p(\boldsymbol{Y} \mid \boldsymbol{X}, \boldsymbol{V}, \boldsymbol{B}, \boldsymbol{\Sigma}).
#' }{
#' L(B, \Sigma | Y, X, V) = p(Y | X, V, B, \Sigma).
#' }
#'
#' When \code{npred > 0}, define the variables \code{Y_star = rbind(Y, y)}, \code{X_star = rbind(X, x)}, and \code{V_star = rbind(cbind(V, w), cbind(t(w), v))}.  Then \code{lmn.suff} calculates summary statistics required to estimate the conditional distribution
#' \deqn{
#' p(\boldsymbol{y} \mid \boldsymbol{Y}, \boldsymbol{X}_\star, \boldsymbol{V}_\star, \boldsymbol{B}, \boldsymbol{\Sigma}).
#' }{
#' p(y | Y, X_star, V_star, B, \Sigma).
#' }
#' The inputs to \code{lmn.suff} in this case are \code{Y = Y}, \code{X = X_star}, and \code{V = V_star}.
#'
#' @example examples/lmn-suff.R
#' @export
lmn.suff <- function(Y, X, V, Vtype, npred = 0) {
  # dimensions of the problem
  n <- nrow(Y)
  q <- ncol(Y)
  N <- n+npred
  # fill out covariate matrix
  if(length(X) == 1) X <- matrix(X, N, 1)
  noBeta <- all(X == 0)
  p <- (!noBeta) * ncol(X)
  # variance type
  Vtype <- .get_Vtype(V, Vtype, n, npred)
  # inner product calculations
  Z <- matrix(0, n, q+p+npred)
  Z[,1:q] <- Y
  if(!noBeta) Z[,q+(1:p)] <- X[1:n,]
  if(Vtype == "full") {
    C <- chol(V[1:n,1:n])
    if(npred > 0) {
      Z[,q+p+(1:npred)] <- V[1:n,n+(1:npred)]
    }
    # IP' * V^{-1} * IP
    IP <- crossprod(backsolve(r = C, x = Z, transpose = TRUE))
    # log|V|
    ldV <- 2 * sum(log(diag(C)))
  } else if(Vtype ==  "diag") {
    IP <- matrix(0, q+p+npred, q+p+npred)
    if(length(V) == 1) {
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)])/V
      ldV <- n * log(V)
    } else {
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)], Z[,1:(q+p)]/V[1:n])
      ldV <- sum(log(V[1:n]))
    }
  } else if(Vtype == "acf") {
    acf <- V
    if(npred > 0) {
      Z[,q+p+(1:npred)] <- toeplitz2(rev(acf[1+1:n]), acf[n+1:npred])
    }
    DL <- .DurbinLevinsonEigen(X = Z, Y = matrix(0),
                               acf = acf[1:n], calcMode = 1)
    IP <- DL$IP
    ldV <- DL$ldV
  } else if(Vtype == "Toeplitz") {
    Tz <- V
    IP <- crossprod(Z, solve(Tz, Z))
    ldV <- determinant(Tz, logarithm = TRUE)
  } else {
    stop("Unrecognized variance type.")
  }
  # convert inner products to sufficient statistics
  S <- IP[1:q,1:q,drop=FALSE]
  if(noBeta) {
    Bhat <- NULL
    T <- NULL
  } else {
    T <- IP[q+(1:p),q+(1:p),drop=FALSE]
    Bhat <- solveV(T, IP[q+(1:p),1:q,drop=FALSE])
    S <- S - IP[1:q,q+(1:p)] %*% Bhat
  }
  # predictive distribution
  if(npred > 0) {
    Ap <- IP[q+p+(1:npred),1:q,drop=FALSE]
    if(noBeta) {
      Xp <- NULL
    } else {
      Xp <- X[n+(1:npred),] - IP[q+p+(1:npred),q+(1:p),drop=FALSE]
    }
    if(Vtype == "full") {
      Vp <- V[n+(1:npred),n+(1:npred),drop=FALSE]
    } else if(Vtype == "diag") {
      if(length(V) == 1) {
        Vp <- diag(V, npred)
      } else {
        Vp <- diag(V[n+(1:npred)], npred)
      }
    } else if(Vtype == "acf") {
      Vp <- stats::toeplitz(acf[1:npred])
    } else {
      stop("Unrecognized variance type.")
    }
    Vp <- Vp - IP[q+p+(1:npred),q+p+(1:npred)]
  }
  # output
  ans <- list(Bhat = Bhat, S = S, T = T, ldV = ldV, n = n, p = p, q = q)
  if(npred > 0) {
    ans <- c(ans, list(Ap = Ap, Xp = Xp, Vp = Vp))
  }
  class(ans) <- "lmn_suff"
  ans
}

# helper functions

.get_Vtype <- function(V, Vtype, n, npred) {
  N <- n + npred
  if(is.numeric(V)) {
    if(length(V) == 1) {
      # scalar
      if(!missing(Vtype) && Vtype != "scalar") {
        stop("Incompatible V and Vtype.")
      } else {
        Vtype <- "diag" # or should it be scalar?
      }
    } else if(is.vector(V) && length(V) == N) {
      # diagonal or acf
      if(missing(Vtype)) {
        stop("Vtype cannot be missing for vector inputs.")
      } else if(!Vtype %in% c("diag", "acf")) {
        stop("Incompatible V and Vtype.")
      }
    } else if(is.matrix(V) && all(dim(V) == N)) {
      # full
      if(!missing(Vtype) && Vtype != "full") {
        stop("Incompatible V and Vtype.")
      } else {
        Vtype <- "full"
      }
    }
  } else if(class(V) == "Toeplitz" && nrow(V) == N) {
    if(!missing(Vtype) && Vtype != "acf") {
      stop("Incompatible V and Vtype.")
    } else {
      Vtype <- "Toeplitz"
    }
    if(npred > 0) stop("npred > 0 for V of class 'Toeplitz' not supported.")
  } else {
    if(!is.numeric(V)) {
      stop("V is of unrecognized type.")
    } else {
      stop("V has invalid dimension(s).")
    }
  }
  Vtype
}

## lmn.suff.full <- function(Y, X, V, npred) {
##   # dimensions of the problem
##   n <- nrow(Y)
##   q <- ncol(Y)
##   N <- n+npred
##   noBeta <- all(X == 0)
##   p <- (!noBeta) * ncol(X)
##   # inner product calculations
##   Z <- matrix(0, n, q+p+npred)
##   Z[,1:q] <- Y
##   if(!noBeta) Z[,q+(1:p)] <- X[1:n,]
##   C <- chol(V[1:n,1:n])
##   if(npred > 0) {
##     Z[,q+p+(1:npred)] <- V[1:n,n+(1:npred)]
##   }
##   # IP' * V^{-1} * IP
##   IP <- crossprod(backsolve(r = C, x = Z, transpose = TRUE))
##   # log|V|
##   ldV <- 2 * sum(log(diag(C)))
##   Vp <- V[n+(1:npred),n+(1:npred),drop=FALSE]
##   list(IP = IP, ldV = ldV, Vp = Vp)
## }


# \describe{
#   \item{\code{Bhat}}{The \code{p x q} matrix \code{(t(X)V^{-1}X)^{-1}t(X)V^{-1}Y}, or \code{NULL} if \code{X = 0}.  The \eqn{p \times q}{p x q} matrix \eqn{\hat{\boldsymbol{B}} = (\boldsymbol{X}'\boldsymbol{V}^{-1}\boldsymbol{X})^{-1}\boldsymbol{X}'\boldsymbol{V}^{-1}\boldsymbol{Y}}{B_hat = (X'V^{-1}X)^{-1}X'V^{-1}Y}.}
#   \item{\code{T}}{The \code{p x p} matrix \code{t(X)V^{-1}X}.}
#   \item{\code{S}}{The \code{q x q} matrix \code{t(Y-X Bhat)V^{-1}(Y-X Bhat)}.}
#   \item{\code{ldV}}{The scalar log-determinant of \code{V}.}
#   \item{\code{n}, \code{p}, \code{q}}{The problem dimensions, namely \code{n = nrow(Y)}, \code{p = nrow(Beta)} (or \code{p = 0} if \code{X = 0}), and \code{q = ncol(Y)}.}
# }
# In addition, when \code{npred > 0} and with \code{x}, \code{w}, and \code{v} defined in \strong{Details}:
# \describe{
#   \item{\code{Ap}}{The \code{npred x q} matrix \code{t(w)V^{-1}Y}.}
#   \item{\code{Xp}}{The \code{npred x p} matrix \code{x - t(w)V^{-1}X}.}
#   \item{\code{Vp}}{The scalar \code{v - t(w)V^{-1}w}.}
# }
