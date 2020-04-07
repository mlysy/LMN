## library(LMN)
source("lmn-testfunctions.R")
context("Marginal/Conditional")

test_that("Lik * Prior = Marg * Cond.", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(-1, 0, 1, 3, 4), q = c(1, 3, 5),
                          Vtype = c("scalar", "diag", "acf", "full"),
                          noSigma = c(TRUE, FALSE),
                          noOmega = c(TRUE, FALSE),
                          noPsi = c(TRUE, FALSE),
                          noPrior = c(TRUE, FALSE))
  # remove noBeta && (noSigma || noOmega)
  case.par <- case.par[with(case.par, {p != 0 | (!noSigma & !noOmega)}),]
  # when noPrior always have noOmega = noPsi = TRUE
  case.par <- case.par[with(case.par, {!noPrior | (noOmega & noPsi)}),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
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
    noOmega <- cp$noOmega
    noPsi <- cp$noPsi
    preSuff <- cp$preSuff
    noPrior <- cp$noPrior
    if(noBeta) noOmega <- TRUE
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
    # set up hyperparameters
    if(!noBeta) {
      Lambda <- rMnorm(p,q)
    } else {
      Lambda <- NULL
    }
    if(!noOmega) {
      Omega <- crossprod(rMnorm(p))
    } else {
      Omega <- ifelse(noBeta, NA, 0)
    }
    if(!noPsi) {
      Psi <- crossprod(rMnorm(q))
    } else {
      Psi <- 0
    }
    if(!noSigma) {
      nu <- ifelse(noPrior, 0, rexp(1,rate=1/q) + q)
    } else {
      nu <- NA
    }
    # setup parameters
    if(!noBeta) {
      Beta <- rMnorm(p,q)
    } else {
      Beta <- NULL
    }
    if(!noSigma) {
      Sigma <- crossprod(rMnorm(q))
    } else {
      Sigma <- diag(q)
    }
    # prior
    prior <- list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
    #if(noPrior) {
    #  lpi <- 0
    #} else {
    lpi <- lMNIW(X = Beta, V = Sigma,
                 prior$Lambda, prior$Omega, prior$Psi, prior$nu)
    #}
    # loglik
    if(noBeta) {
      Mu <- matrix(0,n,q)
    } else {
      Mu <- XR %*% Beta
    }
    ll <- lMnorm(Y, Mu, VR, Sigma)
    # logpost & logmarg
    if(Vtype == "acf") {
      Data <- list(Y = Y, X = XX, V = acf, Vtype = Vtype)
    } else {
      Data <- list(Y = Y, X = XX, V = VV, Vtype = Vtype)
    }
    suff <- do.call(lmn_suff, args = Data)
    if(noPrior) {
      prior <- lmn_prior(p = suff$p, q = suff$q, nu = ifelse(noSigma, NA, 0))
    } else {
      prior <- do.call(lmn_prior, c(list(p = suff$p, q = suff$q), prior))
    }
    ## if(noPrior) prior <- lmn_prior(suff, nu = noSigma)
    ## prior2 <- lmn_prior(suff, nu = noSigma)
    ## post2 <- lmn_post(suff, prior = prior2)
    post <- lmn_post(suff, prior = prior)
    lpm <- lmn_marg(suff, prior = prior, post = post)
    lpc <- lMNIW(Beta, Sigma,
                 post$Lambda, post$Omega, post$Psi, post$nu)
    # smallest of relative and absolute error
    mxd <- abs((lpi+ll) - (lpc+lpm))
    mxd <- min(mxd, mxd/abs(lpi+ll))
    if(calc.diff) {
      MaxDiff[ii] <- mxd
    } else {
      expect_equal(mxd, 0, info = paste0("Case ", ii))
    }
  }
})
