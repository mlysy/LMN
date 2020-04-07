#--- test sufficient statistics -------------------------------------------------
## library(LMN)
source("lmn-testfunctions.R")
context("Predictions")

test_that("Prediction statistics are correctly computed.", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                          Vtype = c("scalar", "diag", "acf", "full"),
                          noSigma = c(TRUE, FALSE),
                          npred = c(1, 5))
  ncases <- nrow(case.par)
  n <- 20
  if(calc.diff) {
    MaxDiff <- data.frame(Ap = rep(NA, ncases), Vp = NA, Xp = NA)
  }
  for(ii in 1:ncases) {
    # get parameters
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    Vtype <- as.character(cp$Vtype)
    noSigma <- cp$noSigma
    npred <- cp$npred
    # set up data
    Y <- rMnorm(n,q)
    if(p == 0) {
      XX <- 0
      XR <- matrix(XX, n+npred, 1)
    } else if(p == -1) {
      XX <- rnorm(1)
      XR <- matrix(XX, n+npred, 1)
    } else {
      XX <- rMnorm(n+npred,p)
      XR <- XX
    }
    acf <- exp(-2*(1:(n+npred))/(n+npred))
    diagV <- rexp(n+npred)
    if(Vtype == "scalar") {
      VV <- acf[1]
      VR <- diag(rep(acf[1], n+npred))
    } else if(Vtype == "diag") {
      VV <- diagV
      VR <- diag(diagV)
    } else {
      VR <- toeplitz(acf)
      VV <- VR
    }
    # calculate with lmn_suff
    if(Vtype == "acf") {
      suff <- lmn_suff(Y = Y, X = XX, V = acf, Vtype = Vtype, npred = npred)
    } else {
      suff <- lmn_suff(Y = Y, X = XX, V = VV, Vtype = Vtype, npred = npred)
    }
    # calculate in R
    X <- XR[1:n,,drop = FALSE]
    V <- VR[1:n,1:n,drop = FALSE]
    if(p != 0) {
      T <- crossprod(X, solve(V, X))
      Bhat <- solve(T, crossprod(X, solve(V, Y)))
    } else {
      T <- 0
      Bhat <- rep(0,q)
    }
    S <- crossprod((Y - X%*%Bhat), solve(V, Y - X%*%Bhat))
    ldV <- ldet(V)
    # check predictions
    if(npred > 0) {
      x <- XR[n+1:npred,]
      w <- VR[1:n,n+1:npred]
      v <- VR[n+1:npred,n+1:npred]
      Ap <- crossprod(w, solve(V, Y))
      if(p != 0) {
        Xp <- x - crossprod(w, solve(V, X))
      }
      Vp <- v - crossprod(w, solve(V, w))
      if(p != 0) {
        if(calc.diff) {
          MaxDiff$Xp[ii] <- max(abs(Xp - suff$Xp))
        } else {
          expect_equal(Xp, suff$Xp)
        }
      }
      if(calc.diff) {
        MaxDiff$Ap[ii] <- max(abs(Ap - suff$Ap))
        MaxDiff$Vp[ii] <- max(abs(Vp - suff$Vp))
      } else {
        expect_equal(Ap, suff$Ap)
        expect_equal(Vp, suff$Vp)
      }
    }
  }
})
