#' Convert list of MNIW parameter lists to vectorized format.
#'
#' Converts a list of return values of multiple calls to [lmn_prior()] or [lmn_post()] to a single list of MNIW parameters, which can then serve as vectorized arguments to the functions in \pkg{mniw}.
#'
#' @param x List of `n` MNIW parameter lists.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`Lambda`}{The mean matrices as an array of size `p x p x n`.}
#'   \item{`Omega`}{The between-row precision matrices, as an array of size `p x p x `.}
#'   \item{`Psi`}{The between-column scale matrices, as an array of size `q x q x n`.}
#'   \item{`nu`}{The degrees-of-freedom parameters, as a vector of length `n`.}
#' }
#' @export
list2mniw <- function(x) {
  # problem dimensions
  p <- nrow(x[[1]]$Lambda)
  q <- ncol(x[[1]]$Lambda)
  n <- length(x)
  nLambda <- p*q
  nOmega <- p*p
  nPsi <- q*q
  # unlist to matrix
  x <- matrix(unlist(x, recursive = TRUE),
              nrow = nLambda+nOmega+nPsi+1, ncol = n)
  # relist in correct format
  list(Lambda = array(x[1:nLambda,], dim = c(p,q,n)),
       Omega = array(x[nLambda+1:nOmega,], dim = c(p,p,n)),
       Psi = array(x[nLambda+nOmega+1:nPsi,], dim = c(q,q,n)),
       nu = x[nLambda+nOmega+nPsi+1,])
}


#' Log-determinant of a variance matrix
#'
#' @param V A variance matrix.
#' @return The log-determinant `log(det(V))`.
#' @noRd
ldet <- function(V) {
  ## determinant(V, logarithm = TRUE)$mod[1]
  solveV(V, x = rep(1, nrow(V)), ldV = TRUE)$ldV
}

#' Log of multi-gamma function.
#'
#' @noRd
lmgamma <- function(x, p) {
  p*(p-1)/4 * log(pi) + sum(lgamma(x + (1-1:p)/2))
}

#' Non-symmetric Toeplitz matrix.
#'
#' @noRd
toeplitz2 <- function(col, row, debug = FALSE) {
  # dimensions
  n <- length(col)
  d <- length(row)
  if(col[1] != row[1]) {
    stop("row[1] and col[1] must be the same.")
  }
  CR <- c(rev(row[-1]), col)
  T <- matrix(NA, n, d)
  if(debug) browser()
  for(ii in 1:d) {
    T[,ii] <- CR[(d-ii)+1:n]
  }
  T
}

#' Names of prior
#'
#' @noRd
.PriorNames <- sort(c("Lambda", "Omega", "Psi", "nu"))

#' Default prior specification
#'
#' @noRd
.DefaultPrior <- function(prior, p, q, noSigma) {
  noBeta <- p == 0
  # extract elements
  Lambda <- prior$Lambda
  Omega <- prior$Omega
  nu <- prior$nu
  Psi <- prior$Psi
  # assign defaults
  if(!noBeta) {
    noBeta <- (length(Omega) == 1) && is.na(Omega)
  }
  if(noBeta) {
    Lambda <- 0
    Omega <- NA
  } else {
    if(is.null(Lambda)) Lambda <- matrix(0,p,q)
    if(is.null(Omega) || all(Omega == 0)) Omega <- matrix(0,p,p)
  }
  if(missing(noSigma) || !noSigma) {
    noSigma <- (length(nu) == 1) && is.na(nu)
  }
  if(noSigma) {
    Psi <- 0
    nu <- NA
  } else {
    if(is.null(nu)) nu <- 0
    if(is.null(Psi) || all(Psi == 0)) Psi <- matrix(0,q,q)
  }
  list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
}

#' Solve method for variance matrices.
#'
#' @param V Variance matrix
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than `solve` when `V` is known to be positive-definite.
#'
#' @noRd
solveV <- function(V, x, ldV = FALSE) {
  ## C <- tryCatch(chol(V), error = function(e) NULL)
  ## if(is.null(C)) browser()
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}
