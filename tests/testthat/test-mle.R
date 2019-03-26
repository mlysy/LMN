#--- test MLE -------------------------------------------------------------------
## library(LMN)
require(numDeriv)
source("lmn-testfunctions.R")
context("MLE")

test_that("Gradient = 0 at MLE.", {
  calc.grad <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                          Vtype = c("scalar", "diag", "acf", "full"),
                          noSigma = c(TRUE, FALSE))
  # remove noBeta && noSigma
  case.par <- case.par[with(case.par, {p != 0 | !noSigma}),]
  ncases <- nrow(case.par)
  n <- 20
  if(calc.grad) {
    MaxGrad <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    # get parameters
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    noBeta <- p == 0
    Vtype <- as.character(cp$Vtype)
    noSigma <- cp$noSigma
    # set up data
    Y <- rMnorm(n,q)
    if(p == 0) {
      XX <- 0
    } else if(p == -1) {
      p <- 1
      XX <- rnorm(1)
      XR <- matrix(XX, n, 1)
    } else {
      XX <- rMnorm(n,p)
      XR <- XX
    }
    acf <- exp(-2*(1:n)/n)
    diagV <- rexp(n)
    if(Vtype == "scalar") {
      VV <- acf[1]
      VR <- diag(rep(acf[1], n))
    } else if(Vtype == "diag") {
      VV <- diagV
      VR <- diag(diagV)
    } else {
      VV <- toeplitz(acf)
      VR <- VV
    }
    # calculate with lmn.suff
    if(Vtype == "acf") {
      suff <- lmn.suff(Y = Y, X = XX, V = acf, Vtype = Vtype)
    } else {
      suff <- lmn.suff(Y = Y, X = XX, V = VV, Vtype = Vtype)
    }
    Beta.hat <- suff$Beta.hat
    Sigma.hat <- suff$S/suff$n
    # check that gradient equals 0
    if(!noBeta && !noSigma) {
      logdens <- function(theta) {
        lMnorm(X = Y, RowV = VR,
               Mu = XR %*% matrix(theta[1:(p*q)],p,q),
               ColV = ltri2Sig(theta[p*q + (1:(q*(q+1)/2))]))
      }
      ld.grad <- grad(logdens, x = c(Beta.hat, Sig2ltri(Sigma.hat)))
    } else if(noBeta && !noSigma) {
      logdens <- function(theta) {
        lMnorm(X = Y, RowV = VR,
               Mu = matrix(0,n,q), ColV = ltri2Sig(theta))
      }
      ld.grad <- grad(logdens, x = Sig2ltri(Sigma.hat))
    } else if(!noBeta && noSigma) {
      logdens <- function(theta) {
        lMnorm(X = Y, RowV = VR,
               Mu = XR %*% matrix(theta[1:(p*q)],p,q),
               ColV = diag(q))
      }
      ld.grad <- grad(logdens, x = Beta.hat)
    } else {
      ld.grad <- NA
    }
    if(calc.grad) {
      MaxGrad[ii] <- max(abs(ld.grad))
    } else {
      expect_equal(max(abs(ld.grad)),0,tolerance = 1e-6)
    }
  }
})
