# set the seed if we're on cran
if(!identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(2022)
}


#' Random matrix of iid normals
#' @param n Number or rows
#' @param p Number of columns.  When omitted defaults to \code{p = n}.
#' @return A \code{n x p} matrix of of iid elements each \code{N(0,1)}.
rMnorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

#' log-density of matrix-normal distribution
lMnorm <- function(X, Mu, RowV, ColV, debug = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  if(debug) browser()
  Z <- X-Mu
  Z <- solve(ColV, t(Z)) %*% solve(RowV, Z)
  lp <- n*p*log(2*pi) + sum(diag(Z))
  -.5 * (lp + n * LMN:::ldet(ColV) + p * LMN:::ldet(RowV))
}

#' log-density of Inverse-Wishart distribution
liWish <- function(V, Psi, nu) {
  d <- nrow(V)
  ans <- (nu+d+1) * LMN:::ldet(V)
  if(!all(Psi == 0) && (nu > d-1)) {
    ans <- ans + sum(diag(solve(V, Psi)))
    ans <- ans - nu * LMN:::ldet(Psi)
    ans <- ans + nu*d * log(2) + 2*LMN:::lmgamma(nu/2, d)
  }
  -0.5 * ans
}

#' log-density of MNIW distribution
#'
#' @details Omits the relevant normalizations to accomodate improper priors, i.e \code{Omega = 0}, \code{Psi = 0}, \code{nu <= nrow(V)-1}.  Also accomodates degenerate forms (pure Mnorm or iWish) via \code{Omega = NA} or \code{nu = NA}.
lMNIW <- function(X, V, Lambda, Omega, Psi, nu) {
  ld <- 0
  if(!is.na(nu)) {
    ld <- ld + liWish(V, Psi, nu)
  }
  if(!is.na(Omega) && !all(Omega == 0)) {
    ld <- ld + lMnorm(X, Lambda, solve(Omega), V)
  }
  ld
}

#' Transformation from variance matrix to lower triangular part and back
Sig2ltri <- function(Sigma) {
  Sigma[lower.tri(Sigma,diag=TRUE)]
}
ltri2Sig <- function(ltri) {
  q <- .5 * (-1 + sqrt(1+8*length(ltri)))
  Sigma <- matrix(0,q,q)
  Sigma[lower.tri(Sigma,diag=TRUE)] <- ltri
  Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
  Sigma
}
