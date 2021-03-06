#--- test sufficient statistics -------------------------------------------------
## library(LMN)
source("lmn-testfunctions.R")
context("Sufficient Statistics")

test_that("Sufficient statistics are correctly computed.", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                          Vtype = c("scalar", "diag", "acf", "Toeplitz", "full"),
                          noSigma = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 20
  if(calc.diff) {
    MaxDiff <- data.frame(Bhat = rep(NA, ncases),
                          T = NA, S = NA, ldV = NA)
  }
  for(ii in 1:ncases) {
    # get parameters
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    Vtype <- as.character(cp$Vtype)
    noSigma <- cp$noSigma
    # set up data
    Y <- rMnorm(n,q)
    if(p == 0) {
      XX <- 0
      XR <- matrix(XX, n, 1)
    } else if(p == -1) {
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
      VR <- toeplitz(acf)
      VV <- VR
    }
    # calculate with lmn_suff
    if(Vtype == "acf") {
      suff <- lmn_suff(Y = Y, X = XX, V = acf, Vtype = Vtype)
    } else if(Vtype == "Toeplitz") {
      suff <- lmn_suff(Y = Y, X = XX, V = Toeplitz$new(acf = acf))
    } else if(Vtype == "diag") {
      suff <- lmn_suff(Y = Y, X = XX, V = VV, Vtype = Vtype)
    } else {
      suff <- lmn_suff(Y = Y, X = XX, V = VV)
    }
    # calculate in R
    X <- XR[1:n,,drop = FALSE]
    V <- VR[1:n,1:n,drop = FALSE]
    if(p != 0) {
      T <- crossprod(X, solve(V, X))
      Bhat <- solve(T, crossprod(X, solve(V, Y)))
    } else {
      T <- NULL
      Bhat <- rep(0,q)
    }
    S <- crossprod((Y - X%*%Bhat), solve(V, Y - X%*%Bhat))
    ldV <- ldet(V)
    # check sufficient statistics
    if(p != 0) {
      if(calc.diff) {
        MaxDiff$Bhat[ii] <- max(abs(Bhat - suff$Bhat))
      } else {
        expect_equal(Bhat, suff$Bhat)
      }
    }
    if(calc.diff) {
      MaxDiff$T[ii] <- max(abs(T - suff$T))
      MaxDiff$S[ii] <- max(abs(S - suff$S))
      MaxDiff$ldV[ii] <- max(abs(ldV - suff$ldV))
    } else {
      expect_equal(T, suff$T)
      expect_equal(S, suff$S)
      expect_equal(ldV, suff$ldV)
    }
  }
})
