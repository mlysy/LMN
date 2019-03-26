#' Calculate the sufficient statistics of an LMN model.
#'
#' @param Y An \code{n x q} matrix of responses.
#' @param X An \code{N x p} matrix of covariates, where \code{N = n + npred} (see \strong{Details}). May also be passed as:
#' \itemize{
#'   \item A scalar, in which case the one-column covariate matrix is \code{X = X * matrix(1, N, 1)}.
#'   \item \code{X = 0}, in which case the mean of \code{Y} is known to be zero, i.e., no regression coefficients are estimated.
#' }
#' @param V,Vtype The between-observation variance specification.  Currently the following options are supported:
#' \enumerate{
#'   \item \code{Vtype = "full"}: \code{V} is an \code{N x N} symmetric positive-definite matrix.
#'   \item \code{Vtype = "diag"}: \code{V} is a vector of length \code{N} such that \code{V = diag(V)}.
#'   \item \code{Vtype = "scalar"}: \code{V} is a scalar such that \code{V = V * diag(N)}.
#'   \item \code{Vtype = "acf"}: \code{V} is either a vector of length \code{N}, or an object of class \code{Toeplitz} such that \code{V = toeplitz(V)}.
#' }
#' For \code{V} specified as a matrix or scalar, \code{Vtype} is deduced automatically and need not be specified.
#' @param npred A nonnegative integer.  If positive, calculates sufficient statistics to make predictions for new responses. See \strong{Details}.
#' @return An S3 object of type \code{lmn_suff}, consisting of a list with elements:
#' \itemize{
#'   \item \code{Beta.hat = (t(X)V^{-1}X)^{-1}t(X)V^{-1}Y}, or \code{NULL} if \code{X = 0}.
#'   \item \code{T = t(X)V^{-1}X}.
#'   \item \code{S = t(Y-X*Betahat)V^{-1}(Y-X*Betahat)}.
#'   \item \code{ldV = log(|V|)}.
#'   \item \code{n = nrow(Y)}.
#'   \item \code{p = nrow(Beta)}, or \code{p = 0} if \code{X = 0}.
#'   \item \code{q = ncol(Y)}.
#' }
#' In addition, when \code{npred > 0}:
#' \itemize{
#'   \item \code{Ap = t(w)V^{-1}Y}.
#'   \item \code{Xp = x - t(w)V^{-1}X}.
#'   \item \code{Vp = v - t(w)V^{-1}w}.
#' }
#'
#' @details
#' The multi-response normal linear regression model is defined as
#' \preformatted{
#' Y ~ MNorm(X Beta, V, Sigma),
#' }
#' where MNorm is the Matrix-Normal distribution, i.e., corresponding to the multivariate normal distribution
#' \preformatted{
#' vec(Y) ~ N( vec(X Beta), Sigma %x% V ),
#' }
#' where \code{vec(Y)} is a vector of length \code{n*q} stacking the elements of \code{Y}, and \code{\%x\%} is the Kronecker product.  The function \code{lmn.suff} returns everything required to calculate the (normalized) probability distribution \code{p(Y | Beta, Sigma, X, V)}.  For mathematical details please see \code{vignette("LMN", package = "LMN")}.
#'
#' When \code{npred > 0}, summary statistics are produced to calculate the conditional (predictive) distribution \code{p(y | Y, XX, VV, Beta, Sigma)}, where
#' \preformatted{
#' YY = rbind(Y, y) ~ MNorm(XX*Beta, VV, Sigma),
#' }
#' \code{VV = rbind(cbind(V, w), cbind(t(w), v))} and \code{XX = rbind(X, x)}.  In this case the inputs to \code{lmn.suff} are \code{Y}, \code{V = VV} and \code{X = XX}. Please see vignette for mathematical details.
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
    Beta.hat <- NULL
    T <- NULL
  } else {
    T <- IP[q+(1:p),q+(1:p),drop=FALSE]
    Beta.hat <- solveV(T, IP[q+(1:p),1:q,drop=FALSE])
    S <- S - IP[1:q,q+(1:p)] %*% Beta.hat
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
  ans <- list(Beta.hat = Beta.hat, S = S, T = T, ldV = ldV, n = n, p = p, q = q)
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
