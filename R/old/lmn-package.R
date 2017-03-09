#-------------------------------------------------------------------------------
#
# preliminary functions for (L)inear (M)odels with (N)uisance parameters
#
# mlysy, dtsar, august 2015
#
#-------------------------------------------------------------------------------

#' Random matrix of iid normals
#'
#' @param n Number or rows
#' @param p Number of columns.  When omitted defaults to \code{p = n}.
#' @return A \code{n x p} matrix of of iid elements each \code{N(0,1)}.
#' @export
rMnorm <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

#' Wrapper for the log-determinant of a matrix
#'
#' @param V A square matrix.
#' @return The log-determinant \code{log(det(V))}.
#' @export
ldet <- function(V) {
  determinant(V, log = TRUE)$mod[1]
}

#' Log of Multi-Gamma Function
#' @export
lmgamma <- function(x, p) {
  p*(p-1)/4 * log(pi) + sum(lgamma(x + (1-1:p)/2))
}

#' Calculate the sufficient statistics of an LMN model.
#'
#' @details
#' The model is defined as
#' Y ~ MNorm(X*Beta, V, Sigma),
#' where MNorm is the Matrix-Normal distribution, i.e.
#' vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).
#' The predictive distribution for new observations y is
#' y | Y ~ y | Y ~ MNorm(Ap + Xp * Beta, Vp, Sigma)
#'
#' @param Y an (n x q) matrix.
#' @param X an (n x p) matrix. can also be passed as
#' \itemize{
#'    \item X = 0: in which case there is no intercept
#'    \item a scalar \neq 0, in which case a scaled intercept X = X * matrix(1, n, 1)
#'    is assumed.
#' }
#' @param V an (n x n) matrix.
#' @param diagV a vector of length n such that V = diag(diagV).  If a scalar, take V = diagV * diag(n).
#' @param acf a vector of length n such that V = toeplitz(acf).
#' @param npred: an integer.  if nonzero, returns sufficient statistics to make
#'    predictions for \code{y | Y, Beta, Sigma}, where
#'    \deqn{\code{YY = rbind(Y, y) ~ MNorm(XX*Beta, VV, Sigma)},}
#'    where \code{VV = rbind(cbind(V, w), cbind(t(w), v))} and \code{XX = rbind(X, x)}.
#'    In this case the inputs are \code{V = VV} and \code{X = XX}.
#' @return a list with the following elements:
#' \itemize{
#'    \item \code{Betahat = (X'V^{-1}X)^{-1}X'V^{-1}Y}
#'    \item \code{T = X'V^{-1}X}
#'    \item \code{S = (Y-X*Betahat)'V^{-1}(Y-X*Betahat)}
#'    \item \code{ldV = log(|V|)}
#'    \item \code{n = nrow(Y)}
#' }
#' and if \code{npred > 0}:
#' \itemize{
#'    \item \code{Ap = w'V^{-1}Y}
#'    \item \code{Xp = x - w'V^{-1}X}
#'    \item \code{Vp = v - w'V^{-1}w}
#' }
#' @export
lmn.suff <- function(Y, X, V, diagV, acf, npred = 0, debug = FALSE) {
  # size of problem
  n <- nrow(Y)
  q <- ncol(Y)
  if(length(X) == 1) {
    X <- matrix(X, n+npred, 1)
  }
  noBeta <- all(X == 0)
  p <- (!noBeta) * ncol(X)
  # inner product calculations
  Z <- matrix(NA, n, q+p+npred)
  Z[,1:q] <- Y
  if(!noBeta) {
    Z[,q+(1:p)] <- X[1:n,]
  }
  if(!missing(V)) {
    if(debug) browser()
    C <- chol(V[1:n,1:n])
    if(npred > 0) {
      Z[,q+p+(1:npred)] <- V[1:n,n+(1:npred)]
    }
    # IP' * V^{-1} * IP
    IP <- crossprod(backsolve(r = C, x = Z, transpose = TRUE))
    # log|V|
    ldV <- 2 * sum(log(diag(C)))
  } else if(!missing(diagV)) {
    if(debug) browser()
    IP <- matrix(0, q+p+npred, q+p+npred)
    if(length(diagV) == 1) {
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)])/diagV
      ldV <- n * log(diagV)
    } else {
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)], Z[,1:(q+p)]/diagV[1:n])
      ldV <- sum(log(diagV[1:n]))
    }
  } else {
    stop("Currently only V and diagV are supported.")
  }
  # convert inner products to sufficient statistics
  S <- IP[1:q,1:q]
  if(noBeta) {
    Beta.hat <- NULL
    T <- 0
  } else {
    T <- IP[q+(1:p),q+(1:p)]
    Beta.hat <- solve(T, IP[q+(1:p),1:q])
    S <- S - IP[1:q,q+(1:p)] %*% Beta.hat
  }
  # predictive distribution
  if(npred > 0) {
    Ap <- IP[q+p+(1:npred),1:q]
    if(noBeta) {
      Xp <- NULL
    } else {
      Xp <- X[n+(1:npred),] - IP[q+p+(1:npred),q+(1:p)]
    }
    if(!missing(V)) {
      Vp <- V[n+(1:npred),n+(1:npred)]
    } else if(!missing(diagV)) {
      if(length(diagV) == 1) {
        Vp <- diag(rep(diagV, npred))
      } else {
        Vp <- diag(diagV[n+(1:npred)])
      }
    } else {
      stop("Currently only V and diagV are supported.")
    }
    Vp <- Vp - IP[q+p+(1:npred),q+p+(1:npred)]
  }
  # output
  ans <- list(Beta.hat = Beta.hat, S = S, T = T, ldV = ldV, n = n)
  if(npred > 0) {
    ans <- c(ans, list(Ap = Ap, Xp = Xp, Vp = Vp))
  }
  ans
}

#' Log-likelihood function for LMN models
#'
#' @export
lmn.loglik <- function(Beta, Sigma, suff, Y, X, V, diagV, acf) {
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, diagV = diagV, acf = acf)
  }
  # get sufficient statistics for likelihood evaluation
  n <- suff$n
  S <- suff$S
  q <- nrow(S)
  ldV <- suff$ldV
  Beta.hat <- suff$Beta.hat
  T <- suff$T
  noBeta <- is.null(Beta.hat)
  # log-likelihood calculation
  if(!noBeta) {
    Z <- Beta-Beta.hat
    S <- S + crossprod(Z, T %*% Z)
  }
  -0.5 * (sum(diag(solve(Sigma,S))) + q*ldV + n*ldet(Sigma) + n*q * log(2*pi))
}

#' Profile Log-Likelihood for the LMN Model
#'
#' @details The model is defined as
#' Y ~ MNorm(X*Beta, V, Sigma),
#' where MNorm is the Matrix-Normal distribution, i.e.
#' vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).
#' @param suff result of call to lmn.suff (i.e. avoid calculations here)
#' @param noSigma logical. if true assumes that Sigma = diag(ncol(Y))
#' @return the calculated profile log-likelihood
#' @export
lmn.prof <- function(suff, Y, X, V, diagV, acf, noSigma = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, diagV = diagV, acf = acf)
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

#' Parameters of the Posterior Conditional MNIW Distribution
#'
#' @export
lmn.post <- function(suff, Y, X, V, diagV, acf, prior = NULL) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, diagV = diagV, acf = acf)
  }
  n <- suff$n
  Betahat <- suff$Betahat
  T <- suff$T
  S <- suff$S
  p <- nrow(Betahat)
  noBeta <- p == 0
  q <- nrow(S)
  # get prior
  nu <- prior$nu
  if(is.null(nu)) nu <- 0
  Psi <- prior$Psi
  if(is.null(Psi)) Psi <- 0
  noSigma <- is.na(nu)
  Lambda <- prior$Lambda
  if(is.null(Lambda)) Lambda <- matrix(0,p,q)
  if(!noBeta) {
    Omega <- prior$Omega
    if(is.null(Omega)) Omega <- matrix(0,p,p)
    noOmega <- all(Omega == 0)
  } else {
    Omega <- NA
    noOmega <- TRUE
  }
  # calculate posterior parameters
  if(noBeta) {
    Lambdahat <- Lambda
    Omegahat <- NA
  } else {
    Omegahat <- T+Omega
    Lambdahat <- solve(T%*%Betahat + Omega%*%Lambda)
  }
  if(noSigma) {
    nuhat <- NA
    Psihat <- S
  } else {
    nuhat <- nu + (n-p*noOmega)
    Psihat <- Psi + S
  }
  if(!noBeta) {
    Psihat <- Psihat + crossprod(Betahat, T%*%Betahat) +
      crossprod(Lambda,Omega%*%Lambda) -
        crossprod(Lambdahat,Omegahat%*%Lambdahat)
  }
  list(Lambda = Lambdahat, Omega = Omegahat, Psi = Psihat, nu = nuhat)
}

#' Marginal Log-Posterior for the LMN Model
#'
#' @param suff The sufficient statistics to be combined with the MNIW conjugate prior.
#' @param post The parameters of the conditional MNIW distribution.  If missing will use \code{prior} to calculate.
#' @param prior The parameters of the prior.  These are required for correct normalization.
#' @export
lmn.marg <- function(suff, Y, X, V, diagV, acf, post, prior) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- lmn.suff(Y = Y, X = X, V = V, diagV = diagV, acf = acf)
  }
  n <- suff$n
  Betahat <- suff$Betahat
  S <- suff$S
  ldV <- suff$ldV
  p <- nrow(Betahat)
  noBeta <- p == 0
  q <- nrow(S)
  # posterior parameters
  if(missing(post)) {
    post <- lmn.post(suff = suff, prior = prior)
  }
  Omegahat <- post$Omega
  Psihat <- post$Psi
  nuhat <- post$nu
  noSigma <- is.na(nu)
  # get prior
  nu <- prior$nu
  if(is.null(nu)) nu <- 0
  Psi <- prior$Psi
  if(is.null(Psi)) Psi <- 0
  noSigma <- is.na(nu)
  #Lambda <- prior$Lambda
  #if(is.null(Lambda)) Lambda <- matrix(0,p,q)
  if(!noBeta) {
    Omega <- prior$Omega
    if(is.null(Omega)) Omega <- matrix(0,p,p)
    noOmega <- all(Omega == 0)
  } else {
    Omega <- NA
    noOmega <- TRUE
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
}


#' (L)inear (M)odel for (N)uisance Parameters
#'
#' @description Efficient profile likelihood and marginal posteriors when nuisance parameters are those of linear regression models.
#' @details The methodology applied to the parameters \code{theta, Beta, Sigma} of models of the form
#' \deq{Y | theta ~ MN( X(theta) * Beta, V(theta), Sigma).}
#' A conditionally conjugate prior for (\code{Beta}, \code{Sigma}) is
#' \deq{(Beta, Sigma) | theta ~ MNIW(Lambda, Omega, Psi, nu).}
#' Note that the hyperparameters can depend on \code{theta}, e.g., \code{Lambda(theta)}, \code{Omega(theta)}, etc.  To achieve some of the restricted forms of the model:
#' \itemize{
#' \item \code{Omega = 0}: Flat prior \code{Beta} ~ 1.
#' \item \code{Omega = NA}: No \code{Beta}, i.e., Y ~ MN(0, V, Sigma).
#' \item \code{nu = NA}: No \code{Sigma}, i.e., Y ~ MN(X * Beta, V, I).
#' }
