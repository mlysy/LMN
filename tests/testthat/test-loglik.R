#--- test loglikelihood --------------------------------------------------------

## library(LMN)
source("lmn-testfunctions.R")
context("Loglikelihood")

test_that("Loglikelihood with precomputations is same as long form.", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 5), q = c(1, 3, 5),
                          Vtype = c("scalar", "diag", "acf", "full"),
                          noSigma = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 20
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    # get parameters
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    noBeta <- p == 0
    Vtype <- as.character(cp$Vtype)
    noSigma <- cp$noSigma
    preSuff <- cp$preSuff
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
    # get sufficient statistics
    ## suff <- lmn_suff(Y = Y, X = XX, V = VV, Vtype = Vtype)
    if(Vtype == "acf") {
      suff <- lmn_suff(Y = Y, X = XX, V = acf, Vtype = Vtype)
    } else {
      suff <- lmn_suff(Y = Y, X = XX, V = VV, Vtype = Vtype)
    }
    # full loglikelihood
    if(!noBeta) {
      Beta <- rMnorm(p,q)
      Mu <- XR %*% Beta
    } else {
      Beta <- NULL
      Mu <- matrix(0,n,q)
    }
    if(!noSigma) {
      Sigma <- crossprod(rMnorm(q))
    } else {
      Sigma <- diag(q)
    }
    llR <- lMnorm(X = Y, Mu = Mu, RowV = VR, ColV = Sigma)
    # calculate with LMN
    lls <- lmn_loglik(Beta = Beta, Sigma = Sigma,
                      suff = suff)
    if(calc.diff) {
      MaxDiff[ii] <- abs(llR - lls)
    } else {
      expect_equal(llR, lls)
    }
  }
})
